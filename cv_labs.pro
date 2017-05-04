#-------------------------------------------------
#
# Project created by QtCreator 2017-02-19T17:24:59
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = cv_labs
TEMPLATE = app

CONFIG += c++14

INCLUDEPATH += "/usr/local/include/gsl/"
LIBS += "/usr/local/lib/*.so"

SOURCES += main.cpp\
        mainwindow.cpp \
    kernel.cpp \
    matrix.cpp \
    utils.cpp \
    scalespace.cpp \
    corner_detectors.cpp \
    descriptor.cpp \
    ransac.cpp \
    haugh.cpp

HEADERS  += mainwindow.h \
    kernel.h \
    matrix.h \
    utils.h \
    scalespace.h \
    corner_detectors.h \
    descriptor.h \
    ransac.h \
    haugh.h

FORMS    += mainwindow.ui
