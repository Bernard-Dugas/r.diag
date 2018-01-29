#
# Used to insert the date of the latest current branch commit
# (with a YYYY-MM-DD format) in the two main DIAGTOOLS
# source files just before generating their
# corrresponding executables.
#
# Author: B.Dugas, January 2018
#
Command=${1:-null}
Sortie_CDF=${2:-cdf2ccc2}
if [[ "${Command}" == "driver.f" || "${Command}" == "cdf2ccc.F90" ]]; then
  Last_Commit_Date=`git log -1 --date=iso --pretty=format:%cd | cut -d- -f1-3`
  sed "s/vdate.*= ['].*[']/vdate = \'${Last_Commit_Date}\'/" ${Command} > ${TMPDIR}/aaaa$$
  sed_result=$?
  if [[ ${sed_result} = 0 && "${Command}" == "driver.f" ]]; then
    /bin/rm -f driver.f
    mv ${TMPDIR}/aaaa$$ driver.f
  elif [[ ${sed_result} = 0 && "${Command}" == "cdf2ccc.F90" ]]; then
    mv ${TMPDIR}/aaaa$$ ${Sortie_CDF}.F90
  else
    echo "Error ${sed_result} in sed/command-change-date"
    exit ${sed_result}
  fi
else
  echo -e "\n Specify either driver.f or cdf2ccc.F90 \n"
  exit 1
fi
exit 0
