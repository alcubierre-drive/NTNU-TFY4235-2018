all:
	@cd Problem1.dir/; make
	@mkdir -p Data/
	@python Problem1.py
	@make todo
	@cd Report; make
todo:
	@grep -n --color TODO Problem1.py
show:
	@cd Report; make show
