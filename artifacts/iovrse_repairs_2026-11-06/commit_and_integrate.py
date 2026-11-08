#!/usr/bin/env python3
"""Commit the reviewed local stack and integrate master/main without force pushes.

Default: print targets only. --execute prepares every integration before pushing.
On a network failure use --resume REPORT to finish the prepared pushes.
Original master/main worktrees are never checked out or reset.
"""
import argparse
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

ORDER = ['antitrust', 'coordination', 'trade', 'vertical', 'financial',
         'competitiontoolbox', 'antitrustBayes']
MATURE = {'antitrust', 'trade', 'competitiontoolbox'}
EXCLUDED = {'.scratch', '.codex', '.agents', '.Rhistory', '.RData', '.Renviron',
            '__pycache__', '.pytest_cache', '.cache'}


def git(repo, *args, data=None, check=True):
    result = subprocess.run(['git', '--literal-pathspecs', '-C', str(repo), *args], input=data,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if check and result.returncode:
        raise RuntimeError(f"git {args[0]} failed in {repo}:\n" +
                           result.stderr.decode(errors='replace') +
                           result.stdout.decode(errors='replace'))
    return result


def value(repo, *args):
    return git(repo, *args).stdout.decode().strip()


def excluded(path):
    return any(p in EXCLUDED or p.startswith('.Rlib') or p.endswith('.Rcheck')
               for p in Path(path).parts)


def validate_origin(repo, name):
    origin = value(repo, 'remote', 'get-url', '--push', 'origin')
    if origin not in (f'git@github.com:luciu5/{name}.git',
                      f'https://github.com/luciu5/{name}.git',
                      f'https://github.com/luciu5/{name}'):
        raise RuntimeError(f'Unexpected origin for {name}; inspect it before publishing.')
    if len(value(repo, 'remote', 'get-url', '--push', '--all', 'origin').splitlines()) != 1:
        raise RuntimeError(f'Multiple push URLs for {name}; refusing ambiguous publication.')
    return origin


def remote_head(repo, branch):
    result = value(repo, 'ls-remote', '--heads', 'origin', f'refs/heads/{branch}')
    return result.split()[0] if result else None


def refresh_pins(repo, resolved):
    desc = repo / 'DESCRIPTION'
    text = desc.read_text()
    pattern = r'^Remotes:[^\n]*(?:\n[ \t]+[^\n]*)*'
    def replace_field(match):
        def replace_pin(pin):
            name = pin.group(2)
            if name in ORDER and name not in resolved:
                raise RuntimeError(f'Dependency {name} has not been prepared yet')
            return pin.group(1) + resolved.get(name, pin.group(3))
        return re.sub(r'(luciu5/([^/@,\s]+)@)([0-9a-fA-F]{40})',
                      replace_pin, match.group(0))
    updated = re.sub(pattern, replace_field, text, flags=re.MULTILINE)
    if updated != text:
        desc.write_text(updated)
    path = repo / 'ci/release-dependencies.json'
    if path.exists():
        pins = json.loads(path.read_text())
        for name in pins:
            if name in ORDER:
                if name not in resolved:
                    raise RuntimeError(f'Dependency {name} has not been prepared yet')
                pins[name] = resolved[name]
        path.write_text(json.dumps(pins, indent=2) + '\n')


def require_manual_vertical_release(repo):
    path = repo / '.github/workflows/R-CMD-check.yaml'
    text = path.read_text()
    old = "    if: github.event_name == 'push' && github.ref == 'refs/heads/main'"
    new = ("    if: github.event_name == 'workflow_dispatch' && "
           "github.ref == 'refs/heads/main' && inputs.publish_release == true")
    if old in text:
        if text.count(old) != 1:
            raise RuntimeError('Ambiguous vertical release job; inspect workflow.')
        text = text.replace(old, new)
    elif new not in text:
        raise RuntimeError('Unrecognized vertical publication guard; inspect workflow.')
    if not re.search(r'^      publish_release:', text, re.MULTILINE):
        anchor = '  workflow_dispatch:\n    inputs:\n'
        if text.count(anchor) != 1:
            raise RuntimeError('Unrecognized vertical workflow_dispatch inputs.')
        text = text.replace(anchor, anchor +
            '      publish_release:\n'
            "        description: 'Explicitly publish a GitHub release after checks pass'\n"
            '        type: boolean\n'
            '        default: false\n'
            '        required: false\n', 1)
    path.write_text(text)


def stage_work(repo):
    staged = git(repo, 'diff', '--cached', '--name-only', '-z').stdout.split(b'\0')
    if any(p and excluded(p.decode()) for p in staged):
        raise RuntimeError(f'Excluded local/private paths are already staged in {repo}; review them first.')
    # Updating tracked files must be separate from adding nonignored new files.
    # A combined explicit add also checks ignored ancestors of tracked files.
    # Leave already-staged changes in place and update only index/worktree diffs.
    tracked = git(repo, 'diff', '--name-only', '-z').stdout
    untracked = git(repo, 'ls-files', '--others', '--exclude-standard',
                    *['--exclude='+name for name in sorted(EXCLUDED)],
                    '--exclude=.Rlib*', '--exclude=*.Rcheck', '-z').stdout
    for paths, flags in [(tracked, ('-u',)), (untracked, ())]:
        selected = [p for p in paths.split(b'\0') if p and not excluded(p.decode())]
        if selected:
            git(repo, '-c', 'core.safecrlf=false', 'add', *flags,
                '--pathspec-from-file=-', '--pathspec-file-nul', data=b'\0'.join(selected)+b'\0')


def commit_work(repo):
    stage_work(repo)
    if git(repo, 'diff', '--cached', '--quiet', check=False).returncode:
        git(repo, 'commit', '-m', 'Repair iovrse workflows and refresh coordinated dependency pins')
    return value(repo, 'rev-parse', 'HEAD')


def write_report(path, report):
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent,
                                     prefix='.iovrse-report-', delete=False) as stream:
        json.dump(report, stream, indent=2)
        stream.write('\n');stream.flush();os.fsync(stream.fileno())
    os.replace(stream.name, path)


def production_differences(repo, source, integrated):
    # A regular R/help/test file's executable bit does not change its contents.
    # Compare exact blob IDs, while retaining file-set and file-type checks.
    roots = ('R', 'NAMESPACE', 'inst/stan', 'man', 'tests', 'DESCRIPTION')
    def files(revision):
        result = {}
        for record in git(repo, 'ls-tree', '-r', '-z', revision, '--', *roots).stdout.split(b'\0'):
            if not record: continue
            metadata, path = record.split(b'\t', 1)
            mode, kind, oid = metadata.split()
            category = b'regular' if mode in (b'100644', b'100755') else mode
            result[path] = (category, kind, oid)
        return result
    left, right = files(source), files(integrated)
    return '\n'.join(path.decode() for path in sorted(left.keys() | right.keys())
                     if left.get(path) != right.get(path))


def push_prepared(report_path):
    report = json.loads(report_path.read_text())
    if not report.get('prepared_complete') or [p['name'] for p in report['packages']] != ORDER:
        raise RuntimeError('Refusing an incomplete integration report.')
    for item in report['packages']:
        if item['target'] != ('master' if item['name'] in MATURE else 'main'):
            raise RuntimeError('Unexpected integration target in recovery report.')
        if not all(re.fullmatch('[0-9a-f]{40}', item[k]) for k in ('source_sha','integrated_sha')):
            raise RuntimeError('Recovery report must contain exact commit SHAs.')
    # Check all destinations before the first push, including on resume.
    for item in report['packages']:
        repo = Path(item['repo']); name = item['name']
        if validate_origin(repo, name) != item['origin']:
            raise RuntimeError(f'Origin changed for {name}')
        if value(repo, 'rev-parse', 'refs/heads/refactor') != item['source_sha']:
            raise RuntimeError(f'Local refactor changed for {name}; inspect the prepared plan.')
        for branch, old, new in [('refactor', item['old_refactor'], item['source_sha']),
                                  (item['target'], item['old_target'], item['integrated_sha'])]:
            if remote_head(repo, branch) not in (old, new):
                raise RuntimeError(f'Remote {name}/{branch} changed; no force push will be attempted.')
    for item in report['packages']:
        repo = Path(item['repo'])
        git(repo, 'push', '--atomic', 'origin',
            item['source_sha'] + ':refs/heads/refactor',
            item['integrated_sha'] + ':refs/heads/' + item['target'])
        item['pushed'] = True
        write_report(report_path, report)
        print('PUSHED', item['name'], 'refactor +', item['target'], item['integrated_sha'][:12], flush=True)
    print('All seven integrated. No force pushes, tags, or releases created by this script.')


def prepare(root, report_path, remote_heads=None, prepared_callback=None):
    manifest = json.loads((root/'antitrust/artifacts/coordinated_release_manifest.json').read_text())
    snapshots = json.loads((root/'antitrust/artifacts/iovrse_repairs_2026-11-06/source-manifest.json').read_text())
    report = {'scope': 'Publish refactor and integrate master/main; no package publication', 'packages': []}
    resolved = {}; targets = {}; origins = {}; source_refs = {}
    for name in ORDER:
        repo = root/name; target = 'master' if name in MATURE else 'main'
        origins[name] = validate_origin(repo, name)
        branch = value(repo, 'branch', '--show-current')
        if branch not in ('refactor', 'main') or (branch == 'main' and name in MATURE):
            raise RuntimeError(f'{name} must be on refactor (or main for a newer package).')
        if value(repo, 'diff', '--name-only', '--diff-filter=U'):
            raise RuntimeError(f'Unresolved conflicts in {name}')
        if remote_heads is None:
            git(repo, 'fetch', 'origin')
            targets[name] = remote_head(repo, target)
            source_refs[name] = remote_head(repo, 'refactor')
        else:
            targets[name] = remote_heads[name][target]
            source_refs[name] = remote_heads[name]['refactor']
        if not targets[name]:
            raise RuntimeError(f'Missing remote target {name}/{target}')
        if value(repo, 'rev-parse', 'refs/remotes/origin/'+target) != targets[name]:
            raise RuntimeError(f'Remote {name}/{target} changed during fetch; retry.')
        base = snapshots[name]['original_head']
        git(repo, 'merge-base', '--is-ancestor', base, 'HEAD')
        if source_refs[name]:
            git(repo, 'merge-base', '--is-ancestor', source_refs[name], 'HEAD')
        if name in MATURE:
            git(repo, 'merge-base', '--is-ancestor', manifest['packages'][name]['sha'], targets[name])

    for name in ORDER:
        repo = root/name; target = 'master' if name in MATURE else 'main'
        if value(repo, 'branch', '--show-current') != 'refactor':
            git(repo, 'switch', '-c', 'refactor')
        refresh_pins(repo, resolved)
        if name == 'vertical':
            require_manual_vertical_release(repo)
        source = commit_work(repo)
        baseline = (manifest['packages'][name]['sha'] if name in MATURE
                    else snapshots[name]['original_head'])
        changed = git(repo, 'diff', '--name-only', '-z', snapshots[name]['original_head'], source).stdout
        paths = [p.decode() for p in changed.split(b'\0') if p]
        workspace = report_path.parent/name
        git(repo, 'worktree', 'add', '--detach', str(workspace), targets[name])
        if name not in MATURE and git(repo, 'merge-base', '--is-ancestor', targets[name], source, check=False).returncode == 0:
            git(workspace, 'merge', '--ff-only', source)
        elif paths:
            patch = git(repo, 'diff', '--binary', baseline, source, '--', *paths).stdout
            (report_path.parent/(name+'.patch')).write_bytes(patch)
            if patch:
                git(workspace, 'apply', '--3way', '--index', data=patch)
                if git(workspace, 'diff', '--cached', '--quiet', check=False).returncode:
                    git(workspace, 'commit', '-m', 'Integrate validated iovrse repairs and dependency pins')
        # Preserve remote-only CI/provenance changes; never silently integrate a
        # different economic implementation than the source we just committed.
        differences = production_differences(repo, source, value(workspace, 'rev-parse', 'HEAD'))
        if differences:
            raise RuntimeError(f'{name} integration differs from the reviewed source in:\n{differences}\nInspect {workspace}; nothing has been pushed yet.')
        integrated = value(workspace, 'rev-parse', 'HEAD')
        item = dict(name=name, repo=str(repo), origin=origins[name], target=target,
                    source_sha=source, integrated_sha=integrated,
                    old_target=targets[name], old_refactor=source_refs[name],
                    worktree=str(workspace), pushed=False)
        resolved[name] = prepared_callback(item) if prepared_callback else integrated
        report['packages'].append(item)
        write_report(report_path, report)
        print('PREPARED',name,target,integrated[:12],flush=True)
    report['prepared_complete'] = True
    write_report(report_path, report)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repositories-root', type=Path, default=Path(__file__).resolve().parents[3])
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--execute',action='store_true')
    modes.add_argument('--resume',type=Path)
    args=parser.parse_args()
    if args.resume:
        if not json.loads(args.resume.read_text()).get('prepared_complete'):
            parser.error('Preparation was incomplete; inspect saved worktrees before continuing.')
        push_prepared(args.resume);return
    for name in ORDER:
        print(f'{name}: commit refactor; integrate origin/{"master" if name in MATURE else "main"}')
    if not args.execute:
        print('Run with --execute to commit and push. Scratch/config/library files are excluded.');return
    path=Path(tempfile.mkdtemp(prefix='iovrse-integrate-'))/'report.json'
    print('Recovery report and integration worktrees:',path.parent,flush=True)
    try:
        prepare(args.repositories_root.resolve(),path)
        push_prepared(path)
    except Exception:
        print('STOPPED. Commits/worktrees are preserved at',path.parent,flush=True)
        if path.exists() and json.loads(path.read_text()).get('prepared_complete'):
            print('Resume the prepared pushes with --resume',path,flush=True)
        else:
            print('Preparation was incomplete; no pushes were attempted. Correct the error and rerun --execute, not --resume.',flush=True)
        raise


if __name__=='__main__':
    main()
