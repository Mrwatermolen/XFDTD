#/bin/bash

OS=$(uname -s)

if [ "$OS" = "Darwin" ]; then
  find ./build/examples -type f -perm +111 -print
  find ./build/examples -type f -perm +111 -print | xargs -I {} sh -c '{}'
elif [ "$OS" = "Linux" ]; then
  find ./build/examples -type f -executable -print
  find ./build/example -type f -executable -print | xargs -I {} sh -c '{}'
else
  echo "OS not supported"
fi
