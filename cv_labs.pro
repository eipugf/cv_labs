#-------------------------------------------------
#
# Project created by QtCreator 2017-02-19T17:24:59
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS+= -fopenmp

TARGET = cv_labs
TEMPLATE = app

CONFIG += c++14

SOURCES += main.cpp\
        mainwindow.cpp \
    kernel.cpp \
    matrix.cpp \
    utils.cpp \
    scalespace.cpp \
    corner_detectors.cpp \
    descriptor.cpp \
    kdthree.cpp

HEADERS  += mainwindow.h \
    kernel.h \
    matrix.h \
    utils.h \
    scalespace.h \
    corner_detectors.h \
    descriptor.h \
    kdthree.h

FORMS    += mainwindow.ui
