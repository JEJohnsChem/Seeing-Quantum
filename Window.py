# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 16:40:30 2019

@author: Barri
"""

import tkinter
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg #, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import numpy as np
x=np.linspace(-10,10,1000)
psi0=np.exp(-x**2/2)
dx=x[1]-x[0]
L=x[-1]-x[0]
f=Figure(figsize=(4,4),constrained_layout=True)
f2=Figure(figsize=(4,4), constrained_layout=True)
f3=Figure(figsize=(3,3))
a=f.add_subplot()
a.plot(x,psi0)
a2=f2.add_subplot()
a3=f3.add_subplot()
wavelength_coefficients=np.zeros((nbasis,1))
def HO():
    nbasis=15
    k=int(k_lowest.get())
    K=np.zeros((nbasis,nbasis))
    V=np.zeros((nbasis,nbasis))
    Vx=1/2*x*x
    a.clear()
    a.plot(x,Vx)
    for i in range(nbasis):
        K[i,i]=(i+1)**2*np.pi**2/(2*L*L)
        for j in range(nbasis):
            V[i,j]=np.sum(2/L*np.sin((i+1)*np.pi*(x+10)/L)*Vx*np.sin((j+1)*np.pi*(x+10)/L))*dx
    H=V+K
    E,V=np.linalg.eig(H)
    
    idx=np.argpartition(E,k)
    psis=np.zeros((len(x),nbasis))
    for i in idx[:k]:
        for j in range(nbasis):
            psis[:,i]=psis[:,i]+np.sqrt(2/L)*np.sin((j+1)*np.pi*(x+10)/L)*V[j,i]
        
        a.plot(x,psis[:,i]+E[i])
        a.text(x=10,y=E[i],s='E = {0:.3f}'.format(np.real(E[i])))
    a.axis([-10,10,-.5,11])
    
    print(E[idx[:k]])
    canvas.draw()

def PIB():
    nbasis=200
    kk=int(k_lowest.get())
    K=np.zeros((nbasis,nbasis))
    V=np.zeros((nbasis,nbasis))
    Vx=0*x
    for i in range(len(x)):
#        if x[i]<0:
#            Vx[i]=200
        if x[i]<-5 or x[i]>5:
            Vx[i]=200
    
    a.clear()
    a2.clear()
    a.plot(x,Vx)
    for i in range(nbasis):
        K[i,i]=(i+1)**2*np.pi**2/(2*L*L)
        for j in range(nbasis):
            V[i,j]=np.sum(2/L*np.sin((i+1)*np.pi*(x+10)/L)*Vx*np.sin((j+1)*np.pi*(x+10)/L))*dx
    H=V+K
    E,W=np.linalg.eig(H)
    
    idx=np.argpartition(E,kk)
    psis=np.zeros((len(x),nbasis))
    for i in idx[:kk]:
        if np.abs(E[i]-np.min(E))<0.0001:
            print(W[0:5,i])
        for j in range(nbasis):
            psis[:,i]=psis[:,i]+np.sqrt(2/L)*np.sin((j+1)*np.pi*(x+10)/L)*W[j,i]
            
        a.plot(x,psis[:,i]+E[i])
        a.text(x=10,y=E[i],s='E = {0:.3f}'.format(np.real(E[i])))
    a.axis([-10,10,-.5,np.max(E[idx[:kk]])+1])
    for i in idx[:kk]:
        if E[i]==np.min(E):
#            a2.plot(np.abs(W[:,i])**2,'+')
            
            for j in range(nbasis):
                k=2*np.pi*(j-nbasis//2)/L
                wavelength_coefficients[j]=np.abs(np.sum(1/np.sqrt(L)*np.exp(1j*k*x)*psis[:,i])*dx)**2
                if k==0:
                    q=2  
                else:
                    wavelength=2*np.pi/k
                    a2.plot(wavelength,wavelength_coefficients[j],'+',color='black')
    a.set_xlabel('Distance',fontsize=14)
    canvas.draw()
#    a2.axis([-L/2,L/2,0,0.02])
    
    canvas2.draw()


def delta():
    nbasis=100
    K=np.zeros((nbasis,nbasis))
    V=np.zeros((nbasis,nbasis))
    k=int(k_lowest.get())
    Vx=0*x
    for i in range(len(x)):
        if x[i]>-.05 and x[i]<0.05:
            Vx[i]=-10
    
    a.clear()
    a.plot(x,Vx)
    for i in range(nbasis):
        K[i,i]=(i+1)**2*np.pi**2/(2*L*L)
        for j in range(nbasis):
            V[i,j]=np.sum(2/L*np.sin((i+1)*np.pi*(x+10)/L)*Vx*np.sin((j+1)*np.pi*(x+10)/L))*dx
    H=V+K
    E,V=np.linalg.eig(H)
    
    idx=np.argpartition(E,k)
    psis=np.zeros((len(x),nbasis))
    for i in idx[:k]:
        for j in range(nbasis):
            psis[:,i]=psis[:,i]+np.sqrt(2/L)*np.sin((j+1)*np.pi*(x+10)/L)*V[j,i]
        a.text(x=10,y=E[i],s='E = {0:.3f}'.format(np.real(E[i])))
        a.plot(x,psis[:,i]+E[i])
    #a.axis([-10,10,-10,11])
    print(E[idx[:k]])
    canvas.draw()

def user_defined():
    if user_input.get():
        print(v.get())
        print(v.get())
        Vx_str=user_input.get()
        nbasis=100
        K=np.zeros((nbasis,nbasis),dtype=complex)
        V=np.zeros((nbasis,nbasis), dtype=complex)
        Vx=eval(Vx_str)
    
        
        a.clear()
        a2.clear()
        if not v.get():
            print('Using Clamped')
            a.plot(x,Vx,color='black')
            for i in range(nbasis):
                K[i,i]=(i+1)**2*np.pi**2/(2*L*L)
                for j in range(nbasis):
                    V[i,j]=np.sum(2/L*np.sin((i+1)*np.pi*(x+10)/L)*Vx*np.sin((j+1)*np.pi*(x+10)/L))*dx
        else:
            PBCx=x#(x-np.min(x))/L*(2*np.pi)
            R=L/2/np.pi
#            Vx_str=Vx_str.replace('x','PBCx')
            Vx=eval(Vx_str)
            print(Vx_str)
            a.plot(x,Vx,color='black')
            for i in range(nbasis):
                k=2*pi*(i-nbasis//2)/L
                K[i,i]=k**2/2
                for j in range(nbasis):
                    kp=2*pi*(j-nbasis//2)/L
                    V[i,j]=1/L*np.sum(np.exp(1j*(k-kp)*x)*Vx)*dx
        H=V+K
        print(np.max(V))
        E,W=np.linalg.eig(H)
        kk=int(k_lowest.get())
        idx=np.argpartition(E,kk)
        print('The {} Lowest Energies are {}'.format(kk, E[idx[:kk]]))
        psis=np.zeros((len(x),nbasis),dtype=complex)
        wavelength_coefficients=np.zeros((nbasis,1))
        for i in idx[:kk]:
            if E[i]==np.min(E):
                a3.clear()
                a3.plot(np.linspace(1,nbasis,nbasis),np.real(W[:,i]),np.linspace(1,nbasis,nbasis),np.imag(W[:,i]))
                f3.suptitle('Coefficients basis functions')
                a3.legend(('Real', 'Imag'))
                canvas3.get_tk_widget().grid(row=10,column=0, pady=20)
                canvas3.draw()
            for j in range(nbasis):
                
                if not v.get():
                    
                    psis[:,i]=psis[:,i]+np.sqrt(2/L)*np.sin((j+1)*np.pi*(x+10)/L)*W[j,i]
                    k=(j+1)*pi/L
                    if k==0:
                        wavelength=2*L
                    else:
                        wavelength=2*np.pi/k
                    if E[i]==np.min(E):
                        wavelength_coefficients[j]=np.real(W[j,i]*np.conj(W[j,i]))
                        a2.plot(wavelength,wavelength_coefficients[j],'+',color='black')
                        a2.plot(-wavelength,wavelength_coefficients[j],'+',color='black')
                else:
                    k=2*pi/L*(j-nbasis//2)
                    if k==0:
                        wavelength=0
                    else:
                        wavelength=2*np.pi/k
                    if E[i]==np.min(E):
                        wavelength_coefficients[j]=W[j,i]*np.conj(W[j,i])
                        a2.plot(wavelength,wavelength_coefficients[j],'+',color='black')
                    psis[:,i]=psis[:,i] + 1/np.sqrt(L)*np.exp(-1j*k*x)*W[j,i]
            if not v.get():
                a.text(x=10,y=E[i],s='E = {0:.3f}'.format(np.real(E[i])))
                a.plot(x,np.real(psis[:,i])/np.max(psis[:,i])+E[i])
                a.plot(x,.5*np.imag(psis[:,i])/np.max(np.abs(psis[:,i]))+E[i],'b')
                a.plot([-10,10],[E[i],E[i]],color='silver',linestyle=':')
                a.axis([-10,10,np.min(E)-1, np.max(E[idx[:kk]])+1])
            else:
                a.text(x=10,y=E[i],s='E = {0:.3f}'.format(np.real(E[i])))
                line1,=a.plot(x,.5*psis[:,i]/np.max(np.abs(psis[:,i]))+E[i],'r')
                line2,=a.plot(x,.5*np.imag(psis[:,i])/np.max(np.abs(psis[:,i]))+E[i],'b')
                a.plot([-10,10],[E[i],E[i]],color='silver',linestyle=':')
                a.axis([np.min(x),np.max(x),np.min(E)-1, np.max(E[idx[:kk]])+1])
                
#        print0*x(psis[:,idx[0]])
        a.set_xlabel('Distance',fontsize=16)
        a.set_ylabel('Energy',fontsize=16)
        a2.set_xlabel('Wavelength',fontsize=16)
        a2.set_ylabel('${|c_\lambda|^2}$', fontsize=16)
        
        canvas.draw()
        canvas2.draw()
        
#canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
##toolbar = NavigationToolbar2TkAgg(canvas, self)
##toolbar.update()
#canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
window = tkinter.Tk()

window.option_add('*font','Ariel 14')
v=tkinter.IntVar()
window.title('Intro to Tkinter')
window.geometry("1200x900")
left_frame=tkinter.Frame(window,bg='white')
left_frame.grid(row=0,column=0)
right_frame=tkinter.Frame(window)
right_frame.grid(row=0,column=1)


canvas=FigureCanvasTkAgg(f,left_frame)
canvas.get_tk_widget().grid(row=1)
canvas2=FigureCanvasTkAgg(f2,left_frame)
canvas2.get_tk_widget().grid(row=3)
canvas3=FigureCanvasTkAgg(f3,right_frame)
tkinter.Label(left_frame,text='Wavefunctions',font=('Ariel',18),bg='white').grid(row=0)
tkinter.Label(left_frame,text='Fourier Transform of Ground State', font=('Ariel',18),bg='white').grid(row=2)
button_free_particle=tkinter.Button(right_frame,text='Free Particle',width=20)
button_pib=tkinter.Button(right_frame,text='Particle in a Box',width=20,command=PIB)
button_harmonic_oscillator=tkinter.Button(right_frame,text='Harmonic Oscillator',width=20, command=HO)
button_delta=tkinter.Button(right_frame,text='Delta Function',width=20, command=delta)
button_user=tkinter.Button(right_frame,text='User Defined',width=20, command = user_defined)
tkinter.Label(right_frame,text='Potential Energy Landscape V(x)').grid(row=0,column=0)
entry_default_text=tkinter.StringVar(window, value='Enter a pythonic, numpy function of x')
user_input=tkinter.Entry(right_frame,textvariable=entry_default_text,font=('Times',10),width=38)

radio_button_PBC=tkinter.Radiobutton(right_frame, text='Periodic Boundary Conditions', variable=v,value=1)
radio_button_Clamped=tkinter.Radiobutton(right_frame,text='Clamped Boundary Conditions', variable=v,value=0)
button_free_particle.grid(row=1,column=0)
button_pib.grid(row=2)
button_harmonic_oscillator.grid(row=3)
button_delta.grid(row=4)
button_user.grid(row=5)
user_input.grid(row=6)
radio_button_Clamped.grid(row=7,column=0, sticky=tkinter.W)
radio_button_PBC.grid(row=7,column=1, sticky=tkinter.E)
k_lowest=tkinter.Entry(right_frame,text='# of Lowest Energy States')
k_lowest.insert(0,5)
tkinter.Label(right_frame,text='Enter the number of states to display').grid(row=8)
k_lowest.grid(row=9)
window.mainloop()
