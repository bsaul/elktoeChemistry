{
  description = "Elktoe Chemistry project flake";
  nixConfig = {
    bash-prompt = "elktoe> ";
  };
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-22.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = nixpkgs.legacyPackages.${system};
    in {
      devShells.default =  pkgs.mkShell {
        nativeBuildInputs = [ 

        ];
        buildInputs = [
          # pkgs.cmake
          pkgs.bashInteractive
          pkgs.harfbuzz
          pkgs.harfbuzz.dev
          pkgs.fribidi
          pkgs.fribidi.devdoc
          pkgs.freetype.dev
          pkgs.openssl.dev
          pkgs.libxml2.dev
          pkgs.util-linux 

          pkgs.R
          pkgs.rPackages.renv
          pkgs.rPackages.languageserver
          # pkgs.rPackages.devtools
          # pkgs.rPackages.open
        ];
      }; 

    });
}
