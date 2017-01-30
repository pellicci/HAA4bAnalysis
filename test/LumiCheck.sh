export PATH=$HOME/.local/bin:/afs/cern.ch/cms/lumi/brilconda-1.1.7/bin:$PATH
pip uninstall -y brilws
pip install --install-option="--prefix=$HOME/.local" brilws
pip show brilws
echo "Software installed"
brilcalc lumi -b "STABLE BEAMS" -i $1 --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/normtag_DATACERT.json -u /fb
