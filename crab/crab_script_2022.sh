#this is not meant to be run locally
#
echo Check if TTY
if [ "`tty`" != "not a tty" ]; then
  echo "YOU SHOULD NOT RUN THIS IN INTERACTIVE, IT DELETES YOUR LOCAL FILES"
else

#echo "ENV..................................."
#env 
#echo "VOMS"
#voms-proxy-info -all
echo "CMSSW BASE, python path, pwd"
echo $CMSSW_BASE 
echo $PYTHON_PATH
echo $PWD 

#echo Found Proxy in: $X509_USER_PROXY

YEAR="2022"

python3 crab_script.py $1 $YEAR #First is job number (provided by CRAB), second arg is added by hand since the scriptArgs param never works...
#echo "=================================================================================== ls gives:"
#ls ./

#xecho "=============================================================================================="

fi
