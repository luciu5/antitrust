"""A permanent: local Git fixtures detect lost remote changes, leaked scratch
files, wrong dependency pins, accidental publication, and unsafe resume pushes.
No GitHub requests or real repository mutations occur in these tests.
"""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('integrate',Path(__file__).with_name('commit_and_integrate.py'))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class IntegrationTests(unittest.TestCase):
    def fixture(self):
        tmp=tempfile.TemporaryDirectory(prefix='iovrse-script-test-');self.addCleanup(tmp.cleanup)
        root=Path(tmp.name);releases={};snapshots={}
        deps={'antitrust':[], 'coordination':['antitrust'], 'trade':['antitrust','coordination'],
              'vertical':['antitrust','trade'], 'financial':[],
              'competitiontoolbox':['antitrust','coordination','trade','vertical'],
              'antitrustBayes':['antitrust','coordination','trade']}
        for name in m.ORDER:
            repo=root/name;repo.mkdir();branch='refactor' if name in m.MATURE else 'main'
            m.git(repo,'init','-b',branch)
            m.git(repo,'config','user.name','Fixture');m.git(repo,'config','user.email','fixture@example.invalid')
            m.git(repo,'config','commit.gpgsign','false')
            (repo/'R').mkdir();(repo/'R/model.R').write_text('old\n')
            (repo/'R/unchanged.R').write_text('unchanged economic code\n')
            (repo/'.gitignore').write_text('.Rlib/\n.Rhistory\n')
            (repo/'DESCRIPTION').write_text('Package: '+name+'\nVersion: 1.0.0\n'+
                ('Remotes: '+',\n    '.join('luciu5/'+p+'@'+'0'*40 for p in deps[name])+'\n' if deps[name] else ''))
            if name in ['vertical','financial','antitrustBayes']:
                pins={'vertical':['coordination'],'financial':['antitrust'],'antitrustBayes':['antitrust','coordination','trade']}[name]
                (repo/'ci').mkdir();(repo/'ci/release-dependencies.json').write_text(json.dumps({k:'0'*40 for k in pins}))
            workflow=repo/'.github/workflows/qa.yml';workflow.parent.mkdir(parents=True);workflow.write_text('old config\n')
            if name=='vertical':
                (repo/'.github/workflows/R-CMD-check.yaml').write_text(
                    'on:\n  workflow_dispatch:\n    inputs:\n      release_certification:\n        type: boolean\n'
                    'jobs:\n  check:\n    runs-on: ubuntu-latest\n  release:\n'
                    "    if: github.event_name == 'push' && github.ref == 'refs/heads/main'\n")
            m.git(repo,'add','.');m.git(repo,'commit','-m','Original local baseline')
            base=m.value(repo,'rev-parse','HEAD');snapshots[name]={'original_head':base}
            bare=root/(name+'.git');m.git(root,'init','--bare',str(bare));m.git(repo,'remote','add','origin',str(bare))
            target='master' if name in m.MATURE else 'main'
            if name in m.MATURE:
                release=root/(name+'-release');m.git(repo,'worktree','add','--detach',str(release),base)
                (release/'R/model.R').write_text('released\n');(release/'.github/workflows/qa.yml').write_text('released config\n')
                m.git(release,'add','.');m.git(release,'commit','-m','Frozen release')
                releases[name]={'sha':m.value(release,'rev-parse','HEAD')}
                (release/'.github/workflows/qa.yml').write_text('released config\nremote mainline selector\n')
                (release/'R/unchanged.R').chmod(0o755)
                m.git(release,'update-index','--chmod=+x','R/unchanged.R')
                m.git(release,'add','.');m.git(release,'commit','-m','Remote-only workflow improvement')
                m.git(release,'push','origin','HEAD:refs/heads/'+target)
                m.git(repo,'push','origin','HEAD:refs/heads/refactor')
                workflow.write_text('released config\n')
            else:m.git(repo,'push','origin','HEAD:refs/heads/'+target)
            (repo/'R/model.R').write_text('repaired\n')
            (repo/'.scratch').mkdir();(repo/'.scratch/private.txt').write_text('never publish\n')
            (repo/'.Rlib').mkdir();(repo/'.Rlib/package.txt').write_text('ignored\n')
            (repo/'.Rhistory').write_text('ignored\n')
            (repo/'literal[notes].txt').write_text('literal path, not a glob\n')
            (repo/':(glob)**').write_text('literal path, not pathspec magic\n')
        a=root/'antitrust/artifacts';a.mkdir()
        (a/'coordinated_release_manifest.json').write_text(json.dumps({'packages':releases}))
        (a/'iovrse_repairs_2026-11-06').mkdir()
        (a/'iovrse_repairs_2026-11-06/source-manifest.json').write_text(json.dumps(snapshots))
        # Resume the exact partially staged state from the failed shell command.
        m.git(root/'antitrust','add','R/model.R')
        report=root/'integration/report.json';report.parent.mkdir()
        return root,report

    def test_complete_push_preserves_mainline_and_pins_and_excludes_private_files(self):
        root,report=self.fixture()
        with patch.object(m,'validate_origin',side_effect=lambda r,n:m.value(r,'remote','get-url','origin')):
            m.prepare(root,report);m.push_prepared(report)
        data=json.loads(report.read_text());resolved={p['name']:p['integrated_sha'] for p in data['packages']}
        for p in data['packages']:
            repo=root/p['name'];sha=p['integrated_sha']
            self.assertEqual(m.remote_head(repo,p['target']),sha)
            self.assertEqual(m.value(repo,'show',sha+':R/model.R'),'repaired')
            self.assertEqual(m.value(repo,'show',sha+':literal[notes].txt'),'literal path, not a glob')
            self.assertEqual(m.value(repo,'show',sha+':'+':(glob)**'),'literal path, not pathspec magic')
            self.assertFalse(any(m.excluded(path) for path in m.value(repo,'ls-tree','-r','--name-only',sha).splitlines()))
            self.assertEqual(m.value(repo,'tag','--list'),'')
            if p['name'] in m.MATURE:
                self.assertIn('remote mainline selector',m.value(repo,'show',sha+':.github/workflows/qa.yml'))
            desc=m.value(repo,'show',sha+':DESCRIPTION')
            for dep,dep_sha in resolved.items():
                if 'luciu5/'+dep+'@' in desc:self.assertIn('luciu5/'+dep+'@'+dep_sha,desc)
        workflow=m.value(root/'vertical','show',resolved['vertical']+':.github/workflows/R-CMD-check.yaml')
        self.assertIn('default: false',workflow);self.assertIn('inputs.publish_release == true',workflow)
        # Network retry: already-pushed refs remain safe and unchanged.
        with patch.object(m,'validate_origin',side_effect=lambda r,n:m.value(r,'remote','get-url','origin')):
            m.push_prepared(report)

    def test_incomplete_report_cannot_push(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'report.json';p.write_text('{"packages": []}')
            with self.assertRaisesRegex(RuntimeError,'incomplete'):m.push_prepared(p)

    def test_production_guard_accepts_only_regular_mode_changes_and_rejects_content_changes(self):
        root,report=self.fixture();repo=root/'antitrust'
        m.commit_work(repo)
        old=m.value(repo,'rev-parse','HEAD')
        m.git(repo,'update-index','--chmod=+x','R/unchanged.R')
        m.git(repo,'commit','-m','Mode-only change')
        mode=m.value(repo,'rev-parse','HEAD')
        self.assertEqual(m.production_differences(repo,old,mode),'')
        (repo/'R/unchanged.R').write_text('different economic code\n')
        m.git(repo,'add','R/unchanged.R');m.git(repo,'commit','-m','Real code change')
        self.assertEqual(m.production_differences(repo,mode,m.value(repo,'rev-parse','HEAD')),'R/unchanged.R')

    def test_tracked_files_inside_newly_ignored_directories_are_updated_without_forcing_new_files(self):
        # A: Git rejects a combined add pathspec containing tracked files below
        # ignored parents. This is the real ai/doc + ai/examples failure.
        with tempfile.TemporaryDirectory(prefix='iovrse-ignored-parent-') as d:
            repo=Path(d);m.git(repo,'init','-b','refactor')
            m.git(repo,'config','user.name','Fixture');m.git(repo,'config','user.email','fixture@example.invalid')
            m.git(repo,'config','commit.gpgsign','false')
            for folder in ['ai/doc','ai/examples']:(repo/folder).mkdir(parents=True)
            (repo/'ai/doc/tracked.md').write_text('unchanged history\n')
            (repo/'ai/examples/tracked.R').write_text('old example\n')
            m.git(repo,'add','.');m.git(repo,'commit','-m','Track files before parent ignore')
            (repo/'.gitignore').write_text('ai/*\n')
            (repo/'ai/examples/tracked.R').write_text('updated example\n')
            (repo/'ai/examples/private-new.R').write_text('ignored and must remain untracked\n')
            (repo/'new.R').write_text('new public source\n')
            m.git(repo,'add','new.R')  # Already staged by the prior failed attempt.
            sha=m.commit_work(repo)
            self.assertEqual(m.value(repo,'show',sha+':ai/examples/tracked.R'),'updated example')
            self.assertEqual(m.value(repo,'show',sha+':ai/doc/tracked.md'),'unchanged history')
            self.assertEqual(m.value(repo,'show',sha+':new.R'),'new public source')
            self.assertNotIn('ai/examples/private-new.R',m.value(repo,'ls-tree','-r','--name-only',sha))

    def test_remote_advance_stops_before_any_publication(self):
        root,report=self.fixture()
        with patch.object(m,'validate_origin',side_effect=lambda r,n:m.value(r,'remote','get-url','origin')):
            m.prepare(root,report)
            data=json.loads(report.read_text());last=data['packages'][-1]
            advance=root/'concurrent';m.git(root/'antitrustBayes','worktree','add','--detach',str(advance),last['old_target'])
            m.git(advance,'commit','--allow-empty','-m','Concurrent remote update')
            m.git(advance,'push','origin','HEAD:refs/heads/main')
            with self.assertRaisesRegex(RuntimeError,'Remote antitrustBayes/main changed'):m.push_prepared(report)
            for item in data['packages']:
                self.assertEqual(m.remote_head(root/item['name'],'refactor'),item['old_refactor'])

    def test_local_branch_advance_stops_before_any_publication(self):
        root,report=self.fixture()
        with patch.object(m,'validate_origin',side_effect=lambda r,n:m.value(r,'remote','get-url','origin')):
            m.prepare(root,report)
            data=json.loads(report.read_text())
            m.git(root/'antitrust','commit','--allow-empty','-m','New local work')
            with self.assertRaisesRegex(RuntimeError,'Local refactor changed'):m.push_prepared(report)
            for item in data['packages']:
                self.assertEqual(m.remote_head(root/item['name'],item['target']),item['old_target'])

    def test_partial_push_failure_resumes_without_rewriting_completed_repositories(self):
        root,report=self.fixture();real_git=m.git
        def interrupted(repo,*args,**kwargs):
            if args and args[0]=='push' and Path(repo).name=='coordination':
                raise RuntimeError('Simulated network interruption')
            return real_git(repo,*args,**kwargs)
        with patch.object(m,'validate_origin',side_effect=lambda r,n:m.value(r,'remote','get-url','origin')):
            m.prepare(root,report)
            with patch.object(m,'git',side_effect=interrupted):
                with self.assertRaisesRegex(RuntimeError,'network interruption'):m.push_prepared(report)
            data=json.loads(report.read_text())
            self.assertTrue(data['packages'][0]['pushed'])
            self.assertFalse(any(p['pushed'] for p in data['packages'][1:]))
            first=data['packages'][0]
            self.assertEqual(m.remote_head(root/'antitrust','master'),first['integrated_sha'])
            m.push_prepared(report)
            self.assertTrue(all(p['pushed'] for p in json.loads(report.read_text())['packages']))

    def test_private_path_filter_covers_ignored_and_unignored_scratch(self):
        for path in ['.scratch/run/fit.rds','.Rlib/package/DESCRIPTION','.Rhistory','.codex/cache','a/__pycache__/b.pyc']:
            self.assertTrue(m.excluded(path),path)
        for path in ['.Rbuildignore','R/model.R','ci/release-dependencies.json']:
            self.assertFalse(m.excluded(path),path)


if __name__=='__main__':unittest.main()
