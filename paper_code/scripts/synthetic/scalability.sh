#!/bin/bash
screen -d -m -S bar$1
screen -S bar$1 -X stuff 'source venv/bin/activate'$(echo -ne '\015')
screen -S bar$1 -X stuff "python3 synthetic.py 0 scalability/$1/0.0/"$(echo -ne '\015')
