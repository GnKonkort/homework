QT += widgets
#QMAKE_CXXFLAGS_RELEASE = $$replace(QMAKE_CXXFLAGS_RELEASE,"-O2","-O3")
#QMAKE_CFLAGS += "-fsanitize=address -fno-omit-frame-pointer -g -O0"
HEADERS       = window.h \
    bessel_approximation.h \
    multipl_nodes_aproximation.h
SOURCES       = main.cpp \
                bessel_approximation.cpp \
                multipl_nodes_aproximation.cpp \
                window.cpp
