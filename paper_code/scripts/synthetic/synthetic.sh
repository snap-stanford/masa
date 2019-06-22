#!/bin/bash
screen -d -m -S foo$1
screen -S foo$1 -X stuff 'source venv/bin/activate'$(echo -ne '\015')
screen -S foo$1 -X stuff "python3 synthetic.py 0 ordered_synthetic_allperturbs/$1/ && python3 synthetic.py 1 ordered_synthetic_allperturbs/$1/"$(echo -ne '\015')
