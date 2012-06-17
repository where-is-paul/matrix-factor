import re 
import os
import subprocess
        
##START = re.compile(r'\[Ordinal/Name Pointer\] Table')                                     
##TABLE = re.compile(r'^\s+\[([\s*\d*)\] (\w*)')                               

START = re.compile(r'ordinal\s+hint\s+RVA\sname')
        
def dump_table(dllname):
    return subprocess.Popen(['dumpbin.exe', '-exports', dllname], stdout=subprocess.PIPE).communicate()[0]
        
def generate_def(dllname, deffile):
    dump = dump_table(dllname).split('\n')
    for i,x in enumerate(dump):
        if 'RVA' in x:
            break        
    else:
        raise ValueError("Symbol table not found")                                        
        
    syms = []
    for line in dump[i+2::]:
        if 'Summary' in line:
            break
        l = line.split();
        if len(l) == 4:
            syms.append((int(l[0]), l[3]))
        
    d = open(deffile, 'w') 
    d.write('LIBRARY        %s\n' % os.path.split(dllname)[-1])              
    d.write(';CODE          PRELOAD MOVEABLE DISCARDABLE\n')                              
    d.write(';DATA          PRELOAD SINGLE\n')                                            
    d.write('\nEXPORTS\n')                                                                
    for s in syms:
        d.write('%s  @%s\n' % (s[1], s[0]))                                             
    d.close()

dllname = r'C:\Windows\System32\python26.dll'
defname = 'python26.def'
generate_def(dllname, defname)
