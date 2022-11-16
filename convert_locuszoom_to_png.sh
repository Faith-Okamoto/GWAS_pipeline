find . -type f -name '*.pdf' -print0 |
  while IFS= read -r -d '' file
    do pdftoppm "${file}" "${file%.*}" -png -f 1 -singlefile
  done
