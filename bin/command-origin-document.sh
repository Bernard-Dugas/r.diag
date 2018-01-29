#
# Used to retrieve the date and comments of the last
# commit in which the source of this r.diag ${Command}
# was modified and inserts them in in the $Header line
# of the source file $Document. Called in the
# .txt recipies of the lspgm Makefile.
#
# Author: B.Dugas, May 2017
#
Command=${1:-null}
Document=${2:-null}
unset RDIAG_GIT_ORIGIN
if [[ -s ${Command} && -s ${Document} ]]; then
  date_info=$(stat -c %y ${Document})
  export RDIAG_GIT_ORIGIN=`git log -1 --date=iso --pretty=format:%f%x20%x40%x20%cd%x20\(%cn\) ${Command}`
  if [[ "${RDIAG_GIT_ORIGIN}" != 'Initial-GIT-version @ 2014-11-03 17:31:51 -0500 (Bernard Dugas)' ]]; then
    # Seulement mettre a jour la doc des fichiers qui ont deja ete modifies sous git
    sed "s/[$]Header:.*[$]/\$Header: ${RDIAG_GIT_ORIGIN} \$/" ${Document} > $TMPDIR/aaaa$$
    sed_result=$?
    if [[ ${sed_result} = 0 ]]; then
      mv $TMPDIR/aaaa$$ ${Document}
      touch -d "${date_info}" ${Document}
    else
      echo "Error ${sed_result} in sed/command-origin-document"
    fi
  fi
else
  echo -e "\n Specify a ptn or ptn90 src/lspgm module name \n"
  exit 1
fi
exit 0
