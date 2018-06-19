##############################################################################
#Copyright:	WANGBO
#Date:		May, 2013

#file:		Makefile (for mesh)
#############################################################################

# SHELL = /bin/sh
cCOMP = g++
cFLAGS=  -Wno-deprecated -framework OpenGL -framework GLUT -framework Foundation -O3 -g -o

#WORK_DIR= /home/wang-bo/work_by_c/lib
WORK_DIR= /Users/zhangrui/Documents/lib
#WORK_DIR= /Users/wangbo/Dropbox/lib
DIR_INC= -I/usr/include -I/opt/local/include/GL
LIB_DIR = -L/opt/local/lib
LIB=  -lgsl -lgslcblas -lglew -lgl -lglu -lglut -lm -lx11

#############################################################

FLAGS = -D _3D_SUPPORTING_ $(cFLAGS)
INCLUDE = $(DIR_INC)

SRC = *.c *.h

OBJS	= main

FORCE = .force


$(OBJS):$(SRC) change
	$(cCOMP) $(FLAGS) $(OBJS) $(SRC) $(INCLUDE) $(LIB_DIR) $(LIB)

change:
	touch main.c

