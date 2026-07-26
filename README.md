# atcoder-env

りすりすの競技プログラミング用環境

## C++20

### 必要なもの

| | 用途 | 備考 |
| --- | --- | --- |
| g++ | ビルド | C++20 が通るもの |
| CMake 3.25 以上 | 手元ビルド | |
| [risundle](https://github.com/TwoSquirrels/risundle) v2.1 以上 | 提出用バンドル | tree-shaking 付き C++ バンドラー |
| [ac-library](https://github.com/atcoder/ac-library) | `<atcoder/all>` | `CPLUS_INCLUDE_PATH` に置く |
| boost | 任意 | 無ければ該当機能が自動で無効になる |
| clangd | 任意 | エディタの補完・診断 |

初回だけ risundle にライブラリを登録する。

```sh
risundle library add risu ./libs/risu
risundle library add ac-library ~/ac-library
```

### ふだんの流れ

```sh
cp templates/cpp20.cc main.cc   # 問題ごとに作り直す
emacs -nw main.cc               # 解く
./run/cpp20.sh                  # 手元で実行 (標準入力にサンプルを貼る)
./bundle/cpp20.sh | clip        # 提出用に 1 ファイル化してクリップボードへ
```

`main.cc` は追跡対象外。CMake は `main.cc` が無ければテンプレートから自動で用意する。

### エディタ

#### clangd を使うエディタ全般

`compile_flags.txt` をリポジトリ直下に置く。パスが環境依存なので追跡対象外にしてある。

```
-std=gnu++20
-Ilibs/risu
-I/path/to/ac-library
-DDEBUG
-Wall
-Wextra
```

`-Ilibs/risu` は `compile_flags.txt` からの相対パスとして解決される。ac-library だけは絶対パスで書く必要がある。`-DDEBUG` を入れておくと `dump` などデバッグ側も補完対象になる。

#### CLion

`CMakeLists.txt` をそのまま開けばよい。

## ライセンス

[MIT License](LICENSE)
