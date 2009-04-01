import Tkinter as tk
from math import pi, cos, sin

# This file is a module that allows the graphing (in a Tk window) of the
# relative locations of various points of an XB. Basically what you do at
# this point is to import the module, throw a window on screen with new_root,
# graph an XB with draw_Xb as often as liked, and kill the root when so
# desired with kill_root

# Note, this is a rewrite of a display function used at an earlier time and 
# idiosyncracies may attributable to this



#Define stuff that is useful for all graphics fellows
zero_loc = (40, 240)
px = 10 #scaling for units to pixels
sw = 2 #spring width
sr = 3 #spring radius
#Choose color depending on whether or not the segment is free
col = lambda boo: boo * "Black" + (not boo) * "Blue"
# Some preliminary useful functions
x2px = lambda xval: zero_loc[0]+xval*px # Convert x to pixs
y2px = lambda yval: zero_loc[1]-yval*px # Convert y to pixs
xytopx = lambda coord: (zero_loc[0]+coord[0]*px, zero_loc[1]-coord[1]*px)


class graphXB():
    def __init__(self):
        #Make a new root window and name it
        self.root = tk.Tk()
        self.root.title("XB")
        self.root.resizable(0, 0)
        #Create a frame to hold the canvas we will draw on
        self.frame = tk.Frame(self.root, bd=3)
        self.frame.pack()
        #Create our canvas
        self.canvas = tk.Canvas(self.frame, width=200, height=250, bd=0, 
                                highlightthickness=0)
        self.canvas.pack()
        #Display the whole thing
        self.root.update()
        
    def draw_XB(self, XB, free=(True, True, True, True), energies=(0,0,0,0)):
        """Draw the XB on the canvas for the first time"""
        # Decide if update is more apropos
        if self.canvas.find_withtag('fil') is not (): #if we've called before
            try:
                self.update_XB(XB, energies) #knock on the other door
                return
            except TclError: # This should keep the shell from hanging when
                return       # the window is force quit, but doesn't
        #scale and translate locs
        c_loc = xytopx(XB.conv_loc)
        h_loc = xytopx(XB.head_loc)
        #draw lines
        self.canvas.create_line(x2px(-3), y2px(0),
                                x2px(15), y2px(0), width=4,
                                tag ="fil")
        self.canvas.create_line(x2px(0), y2px(0),
                                c_loc[0],c_loc[1],
                                fill=col(free[1]), width=sw,
                                tag="n")
        self.canvas.create_line(c_loc[0],c_loc[1],
                                h_loc[0],h_loc[1],
                                fill=col(free[3]), width=sw,
                                tag="g")
        self.canvas.create_oval(x2px(0)-sr, y2px(0)-sr,
                                x2px(0)+sr, y2px(0)+sr,
                                fill=col(free[0]), outline=col(free[0]),
                                tag="t")
        self.canvas.create_oval(c_loc[0]-sr, c_loc[1]-sr,
                                c_loc[0]+sr, c_loc[1]+sr,
                                fill=col(free[2]), outline=col(free[2]),
                                tag="c")
        self.canvas.create_rectangle(h_loc[0]-sr, h_loc[1]-sr,
                                     h_loc[0]+sr, h_loc[1]+sr,
                                     tag="head")
        #get and draw status text at top of canvas
        #xb_text = 'XB Class: '+XB.give_me_curr_xb_class()
        #bad_text = 'Badness: '+'%.3f' % XB.give_me_badness(head_xy)
        #xb_id = canvas.create_text(10, 5, anchor=NW, text=xb_text, tag="xb")
        #bad_id = canvas.create_text(10, 20, anchor=NW, text=bad_text, tag="bad")
        self.canvas.update()
        
    def update_XB(self, XB, energies=(0,0,0,0)):
        #scale and translate locs
        c_loc = xytopx(XB.conv_loc)
        h_loc = xytopx(XB.head_loc)
        #draw lines
        self.canvas.coords('n',zero_loc[0],zero_loc[1],c_loc[0],c_loc[1])
        self.canvas.coords('c',c_loc[0]-sr,c_loc[1]-sr,c_loc[0]+sr,c_loc[1]+sr)
        self.canvas.coords('g',c_loc[0],c_loc[1],h_loc[0],h_loc[1])
        self.canvas.coords('head',h_loc[0]-sr,h_loc[1]-sr,h_loc[0]+sr,h_loc[1]+sr)
        # get and draw status text at top of canvas
        #self.xb_text = 'XB Class: '+XB.give_me_curr_xb_class()
        #self.bad_text = 'Badness: '+'%.2f' % XB.give_me_badness(head_xy)
        #self.canvas.itemconfigure("xb", text=xb_text)
        #self.canvas.itemconfigure("bad", text=bad_text)
        self.canvas.update()
        
    def kill_root(self):
        #destroy the root window containing a canvas of interest
        self.canvas.winfo_toplevel().destroy()


#canvas = new_root() #now done from calling environment

#use root.update_idletasks() to redraw root window
#use root.update() to process events such as window closings
