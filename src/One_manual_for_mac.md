## Try installation with the modified *mac* files 
- install GNU grep via Homebrew
 ```bash
brew install grep
```
- Installation **miniconda**: 
```
bash etc/scripts/install_miniconda_mac.sh
```
- Installation **ROOT**:
```
bash etc/scripts/install_root.sh
```

Compile by running 
```
make -f Makefile_Mac
```

Then follow README.md again 

## Some kind of explanatory note 
- Installation via terminal in vs code
1.  Problem with 
```bash
bash etc/scripts/install_miniconda.sh
```
- ERROR: 
```error 
etc/scripts/install_miniconda.sh: line 19: 
source $HOME/.bashrc
/Users/name/.bashrc: No such file or directory 
etc/scripts/install_miniconda.sh: line 30: /Users/ida/miniconda/bin/conda: No such file or directory
etc/scripts/install_miniconda.sh: line 31: /Users/ida/miniconda/bin/conda: No such file or directory
```
- On macOS, the default shell is usually **zsh** (not bash), and `.bashrc` may not exist by default.
- The script tries to `source $HOME/.bashrc`, but if that file doesn't exist, you get this error.
- solved by running:
```bash
source $HOME/.bashrc 2>/dev/null || source $HOME/.zshrc 2>/dev/null
```
- does the following: 
	- Tries to run `.bashrc` file hich sets up environment variables and shell settings for Bash (`$HOME/.bashrc`). 
	- Any error messages are hidden (`2>/dev/null`) 
	- ||:  If the previous command fails
	-  Tries to source the `.zshrc` file instead, which is the equivalent setup file for Zsh. Again, errors are hidden (`source $HOME/.zshrc 2>/dev/null`).

Then run the following commands from `install_miniconda.sh`,
but replace minconda by miniforge3 in line 30 and 31

2. Problem with executing 
```bash
make 
```
ERROR:
```error
grep: invalid option -- P
usage: grep [-abcdDEFGHhIiJLlMmnOopqRSsUVvwXxZz] [-A num] [-B num] [-C[num]]
        [-e pattern] [-f file] [--binary-files=value] [--color=when]
        [--context[=num]] [--directories=action] [--label] [--line-buffered]
        [--null] [pattern] [file ...]
g++ -I /Users/ida/miniforge3/include/ -I / -fPIC -g -O2 -Wall -Wextra -Wshadow -pedantic -fopenmp -std=c++20 -c src/PMT.cc -o src/PMT.o          
clang++: error: unsupported option '-fopenmp'
make: *** [src/PMT.o] Error 1
```
That means: 
A.  grep -oP (Perl regex) is not supported by macOS's default grep 
A. SOLUTION:
- install GNU grep via Homebrew
 ```bash
brew install grep
```
- Change grep to ggrep in the Makefile line 12
```c++
CXX_STD := $(shell echo $(ROOT_CFLAGS) | ggrep -oP '(?<=-std=c\+\+)\d+')
```
and B. macOS's default clang++ does not support OpenMP - I tryed different options to solve this but here is the one that worked
B. SOLUTION:
- Install libomp (might already be installed)
```bash
brew install libomp
```
- Replace `-fopenmp` in your `CXXFLAGS` with `-Xpreprocessor -fopenmp`, and add `-lomp` to your linker flags.
```c++
# ...existing code...
#line 7 add -Xpreprocessor 
CXXFLAGS        =-fPIC -g -O2 -Wall -Wextra -Wshadow -pedantic -Xpreprocessor -fopenmp
#add afterwards
LDFLAGS         =-lomp
# line 39 add ${LDFLAGS}
all:    ${OBJS} ${DICTO}
        ${LD} -shared ${CXXFLAGS} ${LDFLAGS} -o ${SLL} ${OBJS} ${LIBSLIN}	
        $(MAKE) clean-intermediate
        @echo "Success!"
```
Now I was able to execute make and run the example files