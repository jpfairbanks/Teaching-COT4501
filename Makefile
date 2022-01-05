syllabus.pdf: Syllabus/syllabus.tex
	pdflatex $^

build/Assessment/quiz1.pdf: Assessment/quiz1.md
	pandoc -o $@ $^
build/Assessment/quiz1_makeup.pdf: Assessment/quiz1_makeup.md
	pandoc -o $@ $^
build/Assessment/quiz2.pdf: Assessment/quiz2.md
	pandoc -o $@ $^
build/Assessment/hw1.pdf: Assessment/hw1.md
	pandoc -o $@ $^
build/Assessment/hw2.pdf: Assessment/hw2.md
	pandoc -o $@ $^
build/Assessment/exam1.pdf: Assessment/exam1.md
	pandoc -o $@ $^
build/install.html: install.md
	pandoc -o $@ $^

build/Assessment: build/Assessment/hw1.pdf build/Assessment/hw2.pdf build/Assessment/quiz1.pdf build/Assessment/quiz2.pdf build/Assessment/exam1.pdf

all: syllabus.pdf hw1.pdf build/Assessment