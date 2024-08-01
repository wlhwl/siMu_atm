### How to install

- Make sure you have root installed
- spdlog
  ```
  cd include
  git clone https://github.com/gabime/spdlog.git
  cd spdlog && mkdir build && cd build
  cmake -DCMAKE_INSTALL_PREFIX="$HOME/.local" ..
  make -j && make install
  ```
- then the usual thing
  ```
  mkdir build
  cd build
  cmake ..
  make
  ./app/analysis /Users/meihualin/Projects/trident/noise/noise_ana/data/muon.root
  ```
- make some simple plot
  ```
  cd scripts
  python3 ovelay.py
  ```

### How to make class from  
Hits->MakeClass("TridentHits")
