from pathlib import Path
import argparse,subprocess,time,json,hashlib

root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description='Run named DT mesh regression cases without overwriting inputs.')
parser.add_argument('--case',nargs='+',required=True)
parser.add_argument('--exe',default='build/Release/dt.exe')
parser.add_argument('--output-dir',required=True)
parser.add_argument('--infolevel',type=int,default=1)
parser.add_argument('--allow-changed',action='store_true',help='Measure intentional algorithm changes without enforcing baseline hashes.')
args=parser.parse_args()
cases=json.loads((root/'test/regression_cases.json').read_text(encoding='utf8'))
out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
exe=(root/args.exe).resolve();runs=[]
for name in args.case:
 case=cases[name];mesh=out/(name+'.vtk');cmd=[str(exe),'--input',str(root/case['input'])]
 if 'adpin' in case:cmd+=['--adpin',str(root/case['adpin'])]
 if 'sizecontrol_profile' in case:cmd+=['--sizecontrol-profile',str(root/case['sizecontrol_profile'])]
 cmd+=case['args']+['--infolevel',str(args.infolevel),'--output',str(mesh)]
 start=time.perf_counter()
 with (out/(name+'.log')).open('w',encoding='utf8') as f:p=subprocess.run(cmd,cwd=root,stdout=f,stderr=subprocess.STDOUT,timeout=3600)
 row=dict(case=name,command=cmd,seconds=time.perf_counter()-start,returncode=p.returncode,exe_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
 if 'sizecontrol_profile' in case:row['profile_sha256']=hashlib.sha256((root/case['sizecontrol_profile']).read_bytes()).hexdigest()
 if mesh.exists():row['output_sha256']=hashlib.sha256(mesh.read_bytes()).hexdigest()
 if 'expected_sha256' in case:row['matches_expected']=row.get('output_sha256')==case['expected_sha256']
 if 'reference' in case:row['matches_reference']=row.get('output_sha256')==hashlib.sha256((root/case['reference']).read_bytes()).hexdigest()
 runs.append(row);(out/'regression_runs.json').write_text(json.dumps(runs,indent=2),encoding='utf8')
 print(json.dumps(row),flush=True)
 assert p.returncode==0 and mesh.exists(),row
 assert args.allow_changed or row.get('matches_expected',True),row
