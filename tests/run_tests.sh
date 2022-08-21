#!/bin/bash

args=$@

usage() {
    echo "Usage: $0 <fichier>"
    echo "run the scripts in scripts_test. e.g.:$0 test2 test3"
    exit
}

main() {
    for file in $args
    do
        if [ -f ./scripts_test/$file.sh ]; then
            ./scripts_test/$file.sh
        else
            echo "./scripts_test/$file.sh do not exist"
        fi
    done
}

if [ $# -le 0 ]; then
    usage
else
    main
fi

