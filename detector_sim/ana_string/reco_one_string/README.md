### How to install

- spdlog
  ```
  cd include
  git clone https://github.com/gabime/spdlog.git
  cd spdlog && mkdir build && cd build
  cmake -DCMAKE_INSTALL_PREFIX="$HOME/.local" ..
  make -j && make install
  ```

### How to make class from  
Hits->MakeClass("TridentHits")
