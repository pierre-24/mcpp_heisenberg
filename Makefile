DIRS = ./lib ./include ./tests ./app

lint:
	 cpplint --quiet --recursive $(DIRS)

todo:
	grep -r "TODO(" $(DIRS)