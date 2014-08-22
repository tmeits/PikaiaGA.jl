#/bin/bash

# Simple bash script to clone remote git repositories
cd ~
if [ -d ./julia-vim/ ]; then
echo 'Remove the git clone'
rm -rf ./julia-vim
fi
git clone git://github.com/JuliaLang/julia-vim.git
cd julia-vim
mkdir -p ~/.vim
echo 'Copying in ~/.vim'
cp -R * ~/.vim
echo 'Done!'
