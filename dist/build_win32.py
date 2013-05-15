import os
import sys
from zipfile import ZipFile

sys.path.insert(0,"..\src");
from main import PROG_VERSION

zip_name = "snipper_win32_%s.zip" % PROG_VERSION;
zip_file = ZipFile(zip_name,'w',allowZip64=True);

bad_dirs = [
  '..\\test',
  '..\\bin-win32',
  '..\\dist',
  '..\\.git',
  '..\\pyenv',
  '..\\freeze'
];

bad_exts = [
  '.pyo',
  '.pyc'
];

bad_files = [
  '.project',
  '.pydevproject'           
];

for dirpath, dirnames, filenames in os.walk(".."):
  if True in map(lambda x: x in dirpath,bad_dirs):
    continue;
  
  zip_dirpath = dirpath.replace("..","snipper");
  for f in filenames:
    if f[-4:] in ('.pyo','.pyc'):
      continue;
    
    if f in bad_files:
      continue;
    
    filepath = os.path.join(dirpath,f);
    zippath = os.path.join(zip_dirpath,f);
    
    print "adding %s" % filepath;
    zip_file.write(filepath,zippath);