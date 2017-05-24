Command=${1:-null}
Document=${2:-null}
unset RDIAG_GIT_ORIGIN
if [[ -s ${Command} && -s ${Document} ]]; then
  date_info=$(stat -c %y${Document})
  export RDIAG_GIT_ORIGIN=`git log -1 --pretty=format:%f%x20%x40%x20%aD%x20\(%an\) ${Command}`
  if [[ "${RDIAG_GIT_ORIGIN}" != 'Initial-GIT-version @ Mon, 3 Nov 2014 17:31:51 -0500 (Bernard Dugas)' ]]; then
    # Seulement mettre a jour la doc des fichiers qui ont deja ete modifies sous git
    sed "s/[$]Header:.*[$]/\$Header: ${RDIAG_GIT_ORIGIN} \$/" ${Document} > $TMPDIR/aaaa$$
    mv $TMPDIR/aaaa$$ ${Document}
    touch -d "${date_info}" ${Document}
  fi
else
  echo -e "\n Specify a ptn or ptn90 src/lspgm module name \n"
fi
