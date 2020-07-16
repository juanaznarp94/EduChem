#!/usr/bin/python3
#   Author: Juan Aznar Poveda
#   Technical University of Cartagena, GIT
#   Copyright (C) 2017
# Git repo: https://juanaznarp94@bitbucket.org/juanaznarp94/tfm.git
# -*- coding: utf-8 -*-
# Dual reading

############################################################################################

from var import *
from settings import *
import matplotlib, sys
import time
import math
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
from numpy import arange, sin, pi
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from Tkinter import *
from ttk import *

import Tkinter 
import Image
import ImageTk

cv = False
fv = False
lsv = False

import cmath
import pylab
from decimal import Decimal
import csv
import RPi.GPIO as GPIO
GPIO.setwarnings(False)
GPIO.setmode(GPIO.BOARD)
GPIO.setup(13, GPIO.OUT) 
GPIO.setup(15, GPIO.OUT)

root = Tk()
root.wm_title("eduChem "+ u"\u00A9"+" 2019")
root.geometry('1280x800+100+100')
root.style = Style()

root.style.theme_use('clam')

pb = Progressbar(root, mode='determinate', maximum=300, length=100)
pb.grid(column='1',row='12',columnspan='3',rowspan='1')
pb['length']=460

w = Text(root, width='65', height='12', bg='white', relief = 'groove')
results = Text(root, width='65', height='2', bg='white', relief = 'groove')
#w.grid(column='0',row='12',columnspan='1',rowspan='1')
#results.grid(column='0',row='11',columnspan='1',rowspan='1')
w.grid(column='1',row='10',columnspan='3',rowspan='1')
results.grid(column='1',row='11',columnspan='3',rowspan='1')

w.insert('1.0', 'Welcome to Educational Electrochemical Client \n Interface. Please:\n 1) Insert the Dual SPE in the adapters plug.\n 2) Choose your fit config\n    (TIA and OPMODE).\n 3) Click Start and save graph.\n'+'\n'+'\n'+'Technical University of Cartagena'+'\n'+ 'TIC GIT 2017 '+u"\u00A9")
#results.insert('1.0', 'Determined concentration results\n(mg and M)')

GRAPH_title = Label(root, text='Edit graph title', font='Roboto 12 bold')
k = Entry(root,width='15')
k.grid(column='2',row='5',columnspan='1',rowspan='1')
GRAPH_title.grid(column='2',row='4',columnspan='1',rowspan='1')

interval_title = Label(root, text='Reading Gap, FV', font='Roboto 12 bold')
points = Entry(root,width='15', font='Times 12 bold')
int1 = Entry(root,width='15', font='Roboto 12 bold')
int2 = Entry(root,width='15', font='Roboto 12 bold')
points_label = Label(root, text='Time, FV', font='Roboto 12 bold')
potential_fv_title = Label(root, text='Potential, FV', font='Roboto 12 bold')
potential_fv_entry = Entry(root, width='15', font='Roboto 12 bold')

range_sweep_cv = Label(root, text='End Potential, CV',font='Roboto 12 bold')
end_v_cv = Entry(root,width='15', font='Roboto 12 bold')      
start_v_cv = Entry(root,width='15', font='Roboto 12 bold')

range_sweep_lsv = Label(root, text='End Potential, LSV', font='Roboto 12 bold')
end_v_lsv = Entry(root,width='15', font='Roboto 12 bold')      
start_v_lsv = Entry(root,width='15', font='Roboto 12 bold')


SUBSTANCE_A = Label(root, text='Substance A', font='Roboto 12 bold')
sa = Entry(root,width='15')
sa.grid(column='1',row='9',columnspan='1',rowspan='1')
SUBSTANCE_A.grid(column='1',row='8',columnspan='1',rowspan='1')

SUBSTANCE_B = Label(root, text='Substance B', font='Roboto 12 bold')
sb = Entry(root,width='15')
sb.grid(column='2',row='9',columnspan='1',rowspan='1')
SUBSTANCE_B.grid(column='2',row='8',columnspan='1',rowspan='1')




SWEEPER = ['10110000', '10110010', '10110011', '10110100', '10110101', '10110110', '10110111', '10111000', '10111001', '10111010', '10111011', '10111100', '10111101']

invSWEEPER_dicc = {'10110000':0,
                   '10110010':0.05,
                   '10110011':0.1,
                   '10110100':0.15,
                   '10110101':0.2,
                   '10110110':0.25,
                   '10110111':0.3,
                   '10111000':0.35,
                   '10111001':0.4,
                   '10111010':0.45,
                   '10111011':0.5,
                   '10111100':0.55,
                   '10111101':0.6}

SWEEPER_dicc = {0:'10110000',
                   0.05:'10110010',
                   0.1:'10110011',
                   0.15:'10110100',
                   0.2:'10110101',
                   0.25:'10110110',
                   0.3:'10110111',
                   0.35:'10111000',
                   0.4:'10111001',
                   0.45:'10111010',
                   0.5:'10111011',
                   0.55:'10111100',
                   0.6:'10111101'}


TIA_dicc = {'Default':TIACN_TIAG_DEFAULT_RLOAD_010,
            '2.75 KOhms':TIACN_TIAG_2_75_RLOAD_010,
            '3.5 KOhms':TIACN_TIAG_3_50_RLOAD_010,
            "7 KOhms":TIACN_TIAG_7_00_RLOAD_010,
            "14 KOhms":TIACN_TIAG_14_0_RLOAD_010,
            "35 KOhms":TIACN_TIAG_35_0_RLOAD_010,
            "120 KOhms":TIACN_TIAG_120__RLOAD_010,
            "350 KOhms":TIACN_TIAG_350__RLOAD_010}

TIA_values = {'Default':0,
            '2.75 KOhms':2750,
            '3.5 KOhms':3500,
            "7 KOhms":7000,
            "14 KOhms":14000,
            "35 KOhms":35000,
            "120 KOhms":120000,
            "350 KOhms":350000}

OPMODE_dicc = {"Deep Sleep":MODECN_OP_MODE_DEEPSLEEP,
               "2-Lead GRGC":MODECN_OP_MODE_2LEADGNDC,
               "Standby":MODECN_OP_MODE_STANDBY00,
               "3-Lead AC":MODECN_OP_MODE_3LEADAMPC,
               "Temperature MT-OFF":MODECN_OP_MODE_TEMPMEAOF,
               "Temperature MT-ON":MODECN_OP_MODE_TEMPMEAON}

TITLE = Label(root, text=' EduChem Interface', font='Roboto 24 bold')
INITIALIZE = Label(root, text='INITIALIZE SWEEP')

TIA_label = Label(root,text='Transimpedance gain', font='Roboto 12 bold')
OPMODE_label = Label(root,text='Operation mode', font='Roboto 12 bold')
SUBSTANCE_label = Label(root, text='Substances', font='Roboto 12 bold')
METHOD_label = Label(root, text='Amperometric method', font='Roboto 12 bold')

# IMAGE TRANSPARENCY TREATMENT

def color_cmp(orig, comp):
    tolerance = 5   #  You can adjust color tolerance value
    inside = []
    for i in range(len(orig) - 1):
        if (orig[i] - tolerance) <= comp[i] >= (orig[i] + tolerance):
            inside.append(True)
        else:
            inside.append(False)
    if not inside.count(False):
        return 255
    else:
        return 0

def automask(img):
    backcolor = img.getpixel((0, 0))
    alpha_mask = Image.new('1', img.size)
    alpha_mask.load()
    for y in range(img.size[1]):
        for x in range(img.size[0]):
            pix_color = img.getpixel((x, y))
            if color_cmp(pix_color, backcolor):
                alpha_mask.putpixel((x, y), 1)
            else:
                alpha_mask.putpixel((x, y), 0)
    img.putalpha(alpha_mask)
    return img

############################################################################################
##  PLOT CONFIG AND SWEEP MAIN FUNCTION

f = Figure(figsize=(5,4), dpi=120, facecolor='white', frameon=False,tight_layout=True)
a = f.add_subplot(111,title='Put your title here',
                  xlabel='v, V',
                  ylabel='i, '+ u"\u00B5"+'A',autoscale_on=True)
a.grid(True)

##t = SW_BIAS_N + SW_BIAS_P + SW_BIAS_P[::-1] + SW_BIAS_N[::-1]
t = []
data_a = []
data_b = []

###########

def smooth(y,box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y,box,mode='same')
    return y_smooth

############

#a.plot(t_fv,smooth(data,4),'darkblue')
dataPlot = FigureCanvasTkAgg(f, master=root)
dataPlot.show()
dataPlot.get_tk_widget().grid(column='5',row='3', columnspan='2', rowspan='10')

def on_change(k):
    f = Figure(figsize=(5,4), dpi=120, facecolor='white', frameon=False,tight_layout=True)
    a = f.add_subplot(111,title=k.get(),
    xlabel='Applied voltage (V)',
    ylabel='Registered current ('+ u"\u00B5"+'A)',autoscale_on=True)

def update(dataA, dataB, method):
    if method=="Linear Sweep Voltammetry":
        a.cla()
        a.grid(True)
        a.set_xlabel('v, V')
        a.set_ylabel('i, '+ u"\u00B5"+'A')
        a.set_title(k.get())
        a.plot(t,dataA,'blue', label=sa.get())
        a.plot(t,dataB,'red', label=sb.get())
        a.legend(prop={'size':8})
        dataPlot.draw()
    if method=="Cyclic Voltammetry":
        a.cla()
        a.grid(True)
        a.set_xlabel('v, V')
        a.set_ylabel('i, '+ u"\u00B5"+'A')
        a.set_title(k.get())
        a.plot(t,dataA,'blue', label=sa.get())
        a.plot(t,dataB,'red', label=sb.get())
        a.legend(prop={'size':8})
        dataPlot.draw()
    if method == "Fixed Voltage":
        a.cla()
        a.grid(True)
        a.set_xlabel('t, s')
        a.set_ylabel('i, '+ u"\u00B5"+'A')
        a.set_title(k.get())
        a.plot(t,dataA,'blue', label=sa.get())
        a.plot(t,dataB,'red', label=sb.get())
        a.legend(prop={'size':8})
        a.axvspan(int(int1.get()),int(int1.get())+10, facecolor='#1f77b4', alpha=0.15)
        a.axvspan(int(int2.get()),int(int2.get())+10, facecolor='#1f77b4', alpha=0.15)
        dataPlot.draw()

def printdoutA(num):
    TIAG = TIA_values["{}".format(variable_TIA.get())]
    value = readadcA()
    vref = 2.5
    vmax = 5-(vref/(2**16))
    N = 16
    binmax = ((2**N)-1)
    volts = (vmax*value)/(binmax)+vref
    current = ((volts-(vref/2))/(TIAG))*1000000
    print ">> A/ Step: %5.3f\n >> Voltage: %5.3f V\m >> Current: %5.3f uA" %(num,volts,current)
    return [current, ">> A/ Step: %5.3f\n >> Voltage: %5.3f V\m >> Current: %5.3f uA" %(num,volts,current)]

def printdoutB(num):
    TIAG = TIA_values["{}".format(variable_TIA.get())]
    value = readadcB()
    vref = 2.5
    vmax = 5-(vref/(2**16))
    N = 16
    binmax = ((2**N)-1)
    volts = (vmax*value)/(binmax)+vref
    current = ((volts-(vref/2))/(TIAG))*1000000
    print ">> B/ Step: %5.3f\n >> Voltage: %5.3f V\m >> Current: %5.3f uA" %(num,volts,current)
    return [current, ">> B/ Step: %5.3f\n >> Voltage: %5.3f V\m >> Current: %5.3f uA" %(num,volts,current)]
    
def sweep(TIA,OPMODE):
    global t
    global data_a
    global data_b
    
    GPIO.output(13,True)
    GPIO.output(15,False)
    init(LOCKWR,\
        TIA,\
        REFCN_BIAS_N[0],\
        OPMODE)
    GPIO.output(13,False)
    GPIO.output(15,True)
    init(LOCKWR,\
        TIA,\
        REFCN_BIAS_N[0],\
        OPMODE)
    
    #start_time = time.time()
    #print("--- %s seconds ---" % (time.time() - start_time))
    
    if ("{}".format(variable_METHOD.get())=="Cyclic Voltammetry"):
        GPIO.output(13,False)
        GPIO.output(15,False)
        pb['maximum'] = 2*((float(end_v_cv.get()))/0.05)+1
        pval = 0
        iterator = int(float(end_v_cv.get())/0.05)+2
        UPLIST = SWEEPER[0:iterator]
        DOWNLIST = UPLIST[::-1]
        DOWNLIST = DOWNLIST[1:]
        for p in UPLIST:
            step(p)
            start_time = time.time()
            auxA = printdoutA(invSWEEPER_dicc[p])
            auxB = printdoutB(invSWEEPER_dicc[p])
            t.append(invSWEEPER_dicc[p])
            data_a.append(auxA[0])
            data_b.append(auxB[0])
            stringA = auxA[1]
            stringB = auxB[1]
            pval = pval + 1
            pb['value'] = pval
            pb.update_idletasks()
            w.insert('1.0', stringA +'\n'+'\n')
            w.insert('1.0', stringB +'\n'+'\n')     
            method = "{}".format(variable_METHOD.get())
            update(data_a, data_b, method)
            taux = 1-(time.time() - start_time)
            time.sleep(taux)
            print("--- %s seconds ---" % (time.time() - start_time))
        for p in DOWNLIST:
            step(p)
            start_time = time.time()
            t.append(invSWEEPER_dicc[p])
            auxA = printdoutA(invSWEEPER_dicc[p])
            auxB = printdoutB(invSWEEPER_dicc[p])
            data_a.append(auxA[0])
            data_b.append(auxB[0])
            stringA = auxA[1]
            stringB = auxB[1]
            pval = pval + 1
            pb['value'] = pval
            pb.update_idletasks()
            w.insert('1.0', stringA +'\n'+'\n')
            w.insert('1.0', stringB +'\n'+'\n')     
            method = "{}".format(variable_METHOD.get())
            update(data_a, data_b, method)
            taux = 1-(time.time() - start_time)
            time.sleep(taux)
            print("--- %s seconds ---" % (time.time() - start_time))
        DATAA = data_a
        DATAB = data_b
        
    elif ("{}".format(variable_METHOD.get())=="Fixed Voltage"):
        GPIO.output(13,False)
        GPIO.output(15,False)
        pb['maximum'] = int(points.get())
        pval = 0
        aux = float(potential_fv_entry.get())
        for p in range(int(points.get())):
            step(SWEEPER_dicc[aux])
            start_time = time.time()
            t.append(p)
            auxA = printdoutA(p)
            auxB = printdoutB(p)
            data_a.append(auxA[0])
            data_b.append(auxB[0])
            stringA = auxA[1]
            stringB = auxB[1]
            pval = pval + 1
            pb['value'] = pval
            pb.update_idletasks()
            w.insert('1.0', stringA +'\n'+'\n')
            w.insert('1.0', stringB +'\n'+'\n')
            method = "{}".format(variable_METHOD.get())         
            update(data_a, data_b, method)
            taux = 1-(time.time() - start_time)
            time.sleep(taux)
            print("--- %s seconds ---" % (time.time() - start_time))
        DATAA = data_a
        DATAB = data_b
        
    elif ("{}".format(variable_METHOD.get())=="Linear Sweep Voltammetry"):
        GPIO.output(13,False)
        GPIO.output(15,False)
        pb['maximum'] = ((float(end_v_lsv.get()))/0.05)+1
        pval = 0
        iterator = int(float(end_v_lsv.get())/0.05)+2
        UPLIST = SWEEPER[0:iterator]
        DOWNLIST = UPLIST[::-1]
        DOWNLIST = DOWNLIST[1:]
        for p in UPLIST:
            step(p)
            start_time = time.time()
            auxA = printdoutA(invSWEEPER_dicc[p])
            auxB = printdoutB(invSWEEPER_dicc[p])
            t.append(invSWEEPER_dicc[p])
            data_a.append(auxA[0])
            data_b.append(auxB[0])
            stringA = auxA[1]
            stringB = auxB[1]
            pval = pval + 1
            pb['value'] = pval
            pb.update_idletasks()
            w.insert('1.0', stringA +'\n'+'\n')
            w.insert('1.0', stringB +'\n'+'\n')     
            method = "{}".format(variable_METHOD.get())
            update(data_a, data_b, method)
            taux = 1-(time.time() - start_time)
            time.sleep(taux)
            print("--- %s seconds ---" % (time.time() - start_time))
        DATAA = data_a
        DATAB = data_b
    else:
        w.insert('1.0', 'Choose a right method' +'\n'+'\n')

    pb.stop()

    GPIO.output(13,True)
    GPIO.output(15,False)
    init(LOCKRO,\
        TIACN_TIAG_35_0_RLOAD_010,\
        REFCN_BIAS_N[0],\
        MODECN_OP_MODE_DEEPSLEEP)
    GPIO.output(13,False)
    GPIO.output(15,True)
    init(LOCKRO,\
        TIACN_TIAG_35_0_RLOAD_010,\
        REFCN_BIAS_N[0],\
        MODECN_OP_MODE_DEEPSLEEP)

    return DATAA,DATAB 

############################################################################################
##  TEXT and BUTTONS

def startCV():
    global t
    on_change(k)
    w.insert('1.0', ">> Transimpedance value selectec: {}".format(variable_TIA.get())+'\n'+'\n')
    w.insert('1.0', ">> Operation mode selected: {}".format(variable_OPMODE.get())+'\n'+'\n')
    w.insert('1.0', ">> Starting sweep..."+'\n'+'\n')
    print">> Transimpedance value selectec: {}".format(variable_TIA.get()) 
    print ">> Operation mode selected: {}".format(variable_OPMODE.get())
    print ">> Starting sweep..."
    TIA = TIA_dicc["{}".format(variable_TIA.get())]
    OPMODE = OPMODE_dicc["{}".format(variable_OPMODE.get())]
    SUBSTANCE = "{}".format(variable_SUBSTANCE.get())
    DATAA, DATAB = sweep(TIA,OPMODE)
    method = "{}".format(variable_METHOD.get())
    #update(DATAA, DATAB, method)

    if ("{}".format(variable_METHOD.get())=="Cyclic Voltammetry"):

        strA = "Results/"+k.get()+"_A_CV.csv"
        strB = "Results/"+k.get()+"_B_CV.csv"
        # .CSV CREATION
        csv_outA = open(strA, 'wb')
        mywriterA = csv.writer(csv_outA)
        for row in zip(t, DATAA):
            mywriterA.writerow(row)
        csv_outA.close()
        
        csv_outB = open(strB, 'wb')
        mywriterB = csv.writer(csv_outB)
        for row in zip(t, DATAB):
            mywriterB.writerow(row)
        csv_outB.close()

        w.insert('1.0', "Data exported to .csv succesfully!"+'\n'+'\n')
        
    elif ("{}".format(variable_METHOD.get())=="Fixed Voltage"):

        
        strA = "Results/"+k.get()+"_A_FV.csv"
        strB = "Results/"+k.get()+"_B_FV.csv"

        # AUTOMATIC LEVEL
        mean2A = np.mean(DATAA[int(int2.get()):int(int2.get())+10])
        mean1A = np.mean(DATAA[int(int1.get()):int(int1.get())+10])
        gapA = (mean1A-mean2A)*1000

        mean2B = np.mean(DATAB[int(int2.get()):int(int2.get())+10])
        mean1B = np.mean(DATAB[int(int1.get()):int(int1.get())+10])
        gapB = (mean1B-mean2B)*1000
        
        print ">> Gap P4 (nA): %5.5f" %(gapA)
        w.insert('1.0', ">> Gap P4: %5.5f" %(gapA) +'\n'+'\n')
        print ">> Gap LH (nA): %5.5f" %(gapB)
        w.insert('1.0', ">> Gap LH: %5.5f " %(gapB) +'\n'+'\n')

        # .CSV CREATION
        csv_outA = open(strA, 'wb')
        mywriterA = csv.writer(csv_outA)
        for row in zip(t, DATAA):
            mywriterA.writerow(row)
        csv_outA.close()
        
        csv_outB = open(strB, 'wb')
        mywriterB = csv.writer(csv_outB)
        for row in zip(t, DATAB):
            mywriterB.writerow(row)
        csv_outB.close()
        
        w.insert('1.0', "Data exported to .csv succesfully!"+'\n'+'\n')
        
    elif ("{}".format(variable_METHOD.get())=="Linear Sweep Voltammetry"):

        strA = "Results/"+k.get()+"_A_LSV.csv"
        strB = "Results/"+k.get()+"_B_LSV.csv"
        
        # .CSV CREATION
        csv_outA = open(strA, 'wb')
        mywriterA = csv.writer(csv_outA)
        for row in zip(t, DATAA):
            mywriterA.writerow(row)
        csv_outA.close()
        
        csv_outB = open(strB, 'wb')
        mywriterB = csv.writer(csv_outB)
        for row in zip(t, DATAB):
            mywriterB.writerow(row)
        csv_outB.close()

        w.insert('1.0', "Data exported to .csv succesfully!"+'\n'+'\n')

    #Lock registers
    GPIO.output(13,True)
    GPIO.output(15,False)
    write(1,1)
    GPIO.output(13,False)
    GPIO.output(15,True)
    write(1,1)
    GPIO.output(13,True)
    GPIO.output(15,True)

     #Calibration curves of available substances
    ##############################################
    #Copy this paragraph to add another substance
    #Calibration curve must have logarithmic axis
    if (SUBSTANCE=='Ascorbic Acid'):
            currentA = max(DATAA)
            currentB = max(DATAB)
            ########## DATA TO FILL ############
            W = 176.12 #Molecular weight (g/mol)
            m = 0.9033 #Slope of CC
            b = 4.1644 #Offset of CC
            V = 1
            ####################################
            log_currentA = math.log10(currentA)
            log_MA = (log_currentA-b)/m
            MA = math.pow(10,log_MA)
            mgA = (V*W*MA)/1000
            print ">> A (M): %5.5f\n >> A (mg): %5.5f" %(MA,mgA)
            results.insert('1.0', ">> A (M): %5.5f  \n >> A (mg): %5.5f " %(MA,mgA))
            
            log_currentB = math.log10(currentB)
            log_MB = (log_currentB-b)/m
            MB = math.pow(10,log_MB)
            mgB = (V*W*MB)/1000
            print ">> B (M): %5.5f\n >> B (mg): %5.5f" %(MB,mgB)
            results.insert('1.0', ">> B (M): %5.5f  \n >> B (mg): %5.5f " %(MB,mgB))


    ##############################################
    if (SUBSTANCE=='Other'):
            currentA = max(DATAA)
            currentB = max(DATAB)
            ########## DATA TO FILL ############
            W = 176.12 #Molecular weight (g/mol)
            m = 0.9033 #Slope of CC
            b = 4.1644 #Offset of CC
            V = 1
            ####################################
            log_currentA = math.log10(currentA)
            log_MA = (log_currentA-b)/m
            MA = math.pow(10,log_MA)
            mgA = (V*W*MA)/1000
            print ">> A (M): %5.5f\n >> A (mg): %5.5f" %(MA,mgA)
            results.insert('1.0', ">> A (M): %5.5f  \n >> A (mg): %5.5f " %(MA,mgA))

            log_currentB = math.log10(currentB)
            log_MB = (log_currentB-b)/m
            MB = math.pow(10,log_MB)
            mgB = (V*W*MB)/1000
            print ">> B (M): %5.5f\n >> B (mg): %5.5f" %(MB,mgB)
            results.insert('1.0', ">> B (M): %5.5f  \n >> B (mg): %5.5f " %(MB,mgB))

    ##############################################

def clearCV():
    global t
    global data_a
    global data_b
    a.cla()
    a.grid(True)
    a.set_xlabel('Applied voltage (V)')
    a.set_ylabel('Registered current ('+ u"\u00B5"+'A)')
    a.set_title(k.get())
    data_a = []
    data_b = []
    t = []
    a.plot(t,data_a,'blue')
    dataPlot.show()
    init(LOCKWR,\
        TIACN_TIAG_35_0_RLOAD_010,\
        REFCN_BIAS_N[0],\
        MODECN_OP_MODE_DEEPSLEEP)
    ## w.insert('1.0', "Graph removed succesfully!"+'\n'+'\n')
    print('Graph removed succesfully!')
    return 0

def saveCV():
    str =  "Results/"+k.get()+".png"
    f.savefig(str, dpi=None, facecolor='w', edgecolor='w',
              orientation='landscape', format='png',
              transparent=False, bbox_inches=None, pad_inches=0.1,
              frameon=False)
    w.insert('1.0', "Graph saved succesfully!"+'\n'+'\n')

def exportCV():
    if ("{}".format(variable_METHOD.get())=="Cyclic Voltammetry"):
        str = k.get()+".csv"
        csv_out = open(str, 'wb')
        mywriter = csv.writer(csv_out)
        for row in zip(t_cv, DATA):
            mywriter.writerow(row)
        csv_out.close()
        w.insert('1.0', "Data exported to .csv succesfully!"+'\n'+'\n')
    elif ("{}".format(variable_METHOD.get())=="Fixed Voltage"):
        str = k.get()+".csv"
        csv_out = open(str, 'wb')
        mywriter = csv.writer(csv_out)
        for row in zip(t_fv, DATA):
            mywriter.writerow(row)
        csv_out.close()
        w.insert('1.0', "Data exported to .csv succesfully!"+'\n'+'\n')
        
def closeCV():
    init(LOCKWR,\
        TIACN_TIAG_35_0_RLOAD_010,\
        REFCN_BIAS_N[0],\
        MODECN_OP_MODE_DEEPSLEEP)
    print('Bye!')
    root.destroy()
    return 0

############################################################################################
##  MAIN MENU

menubar = Menu(root)
menubar.add_command(label="Start", command=startCV)
menubar.add_command(label="Clear",command=clearCV)
menubar.add_command(label="Save graph", command=saveCV)
menubar.add_command(label="Save data to DB")
menubar.add_command(label="Close", command=closeCV)

############################################################################################
##  OPTIONS MENU

variable_TIA = StringVar(root)
variable_TIA.set("35 KOhms")
variable_OPMODE = StringVar(root)
variable_OPMODE.set("3-Lead AC")
variable_SUBSTANCE = StringVar(root)
variable_SUBSTANCE.set("None")
variable_METHOD = StringVar(root)
variable_METHOD.set("None")

SUBSTANCE = OptionMenu(root, variable_SUBSTANCE, "None", "Ascorbic Acid")
METHOD = OptionMenu(root, variable_METHOD, "Linear Sweep Voltammetry", "Linear Sweep Voltammetry", "Cyclic Voltammetry", "Fixed Voltage")

TIA = OptionMenu(root, variable_TIA, "35 KOhms", "2.75 KOhms",
                 "3.5 KOhms", "7 KOhms", "14 KOhms", "35 KOhms", "120 KOhms", "350 KOhms")
OPMODE = OptionMenu(root, variable_OPMODE, "3-Lead AC", "Deep Sleep", "2-Lead GRGC", "3-Lead AC",
                 "Standby", "Temperature MT-OFF", "Temperature MT-ON")

def option_changed_SUBSTANCE(*args):
    print "Substance selected: {}".format(variable_SUBSTANCE.get())
    w.insert('1.0', ">> Substance selected: {}".format(variable_SUBSTANCE.get())+'\n'+'\n')

def option_changed_TIA(*args):
    print "Transimpedance value selected: {}".format(variable_TIA.get())
    w.insert('1.0', ">> TIA value selected: {}".format(variable_TIA.get())+'\n'+'\n')
           
def option_changed_OPMODE(*args):
    print "Operation mode selected: {}".format(variable_OPMODE.get())
    w.insert('1.0', ">> Operation mode selected: {}".format(variable_OPMODE.get())+'\n'+'\n')

def option_changed_METHOD(*args):
    global fv
    global cv
    global lsv
    global int1
    global int2
    global points
    global interval_title
    global range_sweep_lsv
    global start_v_lsv
    global end_v_lsv
    global range_sweep_cv
    global start_v_cv
    global end_v_cv
    
    print "Method selected: {}".format(variable_METHOD.get())
    w.insert('1.0', ">> Method selected: {}".format(variable_METHOD.get())+'\n'+'\n')
    if ("{}".format(variable_METHOD.get())=="Fixed Voltage"):
        if (cv == True):
            range_sweep_cv.grid_forget()
            end_v_cv.grid_forget()
        if (lsv == True):
            range_sweep_lsv.grid_forget()
            end_v_lsv.grid_forget()
        interval_title.grid(column='1',row='6',columnspan='2',rowspan='1')
        int1.grid(column= '1',row='7',columnspan='1',rowspan='1')
        int2.grid(column='2',row='7',columnspan='1',rowspan='1')
        interval_title.grid(column='1',row='6',columnspan='2',rowspan='1')
        points_label.grid(column='3',row='6',columnspan='1',rowspan='1') 
        points.grid(column='3',row='7',columnspan='1',rowspan='1')
        potential_fv_title.grid(column='3',row='4',columnspan='1',rowspan='1')
        potential_fv_entry.grid(column='3',row='5',columnspan='1',rowspan='1')

        fv = True
        cv = False
        lsv = False
        
    elif ("{}".format(variable_METHOD.get())=="Cyclic Voltammetry"):
        if (fv == True):
            interval_title.grid_forget()
            int1.grid_forget()
            int2.grid_forget()
            points_label.grid_forget()
            points.grid_forget()
            potential_fv_title.grid_forget()
            potential_fv_entry.grid_forget()
        if (lsv == True):
            range_sweep_lsv.grid_forget()
            end_v_lsv.grid_forget()
        range_sweep_cv.grid(column='3',row='4',columnspan='1',rowspan='1')
        end_v_cv.grid(column='3',row='5',columnspan='1',rowspan='1')
        cv = True
        lsv = False
        fv = False
        
    elif ("{}".format(variable_METHOD.get())=="Linear Sweep Voltammetry"):
        if (fv == True):
            interval_title.grid_forget()
            int1.grid_forget()
            int2.grid_forget()
            points_label.grid_forget()
            points.grid_forget()
            potential_fv_title.grid_forget()
            potential_fv_entry.grid_forget()
        if (cv == True):
            range_sweep_cv.grid_forget()
            end_v_cv.grid_forget()
        range_sweep_lsv.grid(column='3',row='4',columnspan='1',rowspan='1')
        end_v_lsv.grid(column='3',row='5',columnspan='1',rowspan='1')
        lsv = True
        fv = False
        cv = False
                   
variable_TIA.trace("w", option_changed_TIA)
variable_OPMODE.trace("w", option_changed_OPMODE)
variable_SUBSTANCE.trace("w", option_changed_SUBSTANCE)
variable_METHOD.trace("w", option_changed_METHOD)


############################################################################################
##  GRID CONFIG

root.grid_columnconfigure(0, minsize=60)

root.grid_columnconfigure(1, minsize=100)
root.grid_columnconfigure(2, minsize=100)
root.grid_columnconfigure(3, minsize=100)
root.grid_columnconfigure(4, minsize=60)


TITLE.grid(column='1',row='1',columnspan='3',rowspan='1')

SUBSTANCE_label.grid(column='1',row='2',columnspan='1',rowspan='1')
SUBSTANCE.grid(column='1',row='3',columnspan='1',rowspan='1')

METHOD_label.grid(column='2',row='2',columnspan='1',rowspan='1')
METHOD.grid(column='2',row='3',columnspan='1',rowspan='1')

TIA_label.grid(column='3',row='2',columnspan='1',rowspan='1')
TIA.grid(column='3',row='3',columnspan='1',rowspan='1')

OPMODE_label.grid(column='1',row='4',columnspan='1',rowspan='1')
OPMODE.grid(column='1',row='5',columnspan='1',rowspan='1')

LOGO = PhotoImage(file='ait2.png')
LOGO = LOGO.subsample(7,7)


root.configure(background='white')
SUBSTANCE_A.configure(background='white')
SUBSTANCE_B.configure(background='white')
range_sweep_lsv.configure(background='white')
range_sweep_cv.configure(background='white')
points_label.configure(background='white')
int1.configure(background='white')
int2.configure(background='white')
OPMODE_label.configure(background='white')
METHOD_label.configure(background='white')
SUBSTANCE_label.configure(background='white')
TITLE.configure(background='white')
TIA_label.configure(background='white')
GRAPH_title.configure(background='white')
interval_title.configure(background='white') 
points_label.configure(background='white') 
potential_fv_title.configure(background='white')


#masked = automask(Image.open('ait.gif'))
#maskedtk = ImageTk.PhotoImage(masked)

Label(root, image=LOGO).grid(column='5',row='1',columnspan='2',rowspan='1')
#GRAPH.grid(column='2',row='2', columnspan='2', rowspan='7')
root.config(menu=menubar)
root.mainloop()                 


