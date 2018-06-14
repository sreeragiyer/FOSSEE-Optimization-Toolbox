#!/bin/bash


OP="$(git rev-parse HEAD)"
if [ $? != 0 ]; then 
	echo "Could not fetch Git commit id"

else	
	a='#define GIT_ID "'
	b='"'
	c=$a$OP$b

	echo "Git commit id added to fotConfig.h"

	sed -i "/GIT_ID/c\ $c" ./sci_gateway/cpp/fotConfig.h

fi


