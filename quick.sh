function update()
{
	git add quick.sh
    git add .gitignore
    git add clue_performance/plugins/clue_performance.cc
    git add clue_performance/python/task1_cfi.py
    
    git commit -m "Minor update"
    git push origin
}

update
