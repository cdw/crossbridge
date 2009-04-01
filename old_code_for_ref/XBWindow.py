from Tkinter import *
from math import pi, cos, sin

# This file is a module that allows the graphing (in a Tk window) of the
# relative locations of various points of an XB. Basically what you do at
# this point is to import the module, throw a window on screen with new_root,
# graph an XB with draw_Xb as often as liked, and kill the root when so
# desired with kill_root

# TODO:
# * Rewrite draw XB to use an actual XB data structure and pull out values
#   for what is fixed and what have you. DONE
# * Figure out how to call updates with event handlers? DONE (USE .update())
# * What is the diff between root.withdraw() and root.destroy()
# * Add in text display function. DONE
# * Make text display update. DONE



#Define stuff that is useful for all graphics fellows
zero_loc = (40, 240)
px = 10 #scaling for units to pixels
sw = 2 #spring width
sr = 3 #spring radius
#Choose color depending on whether or not the segment is free
col = lambda boo: boo * "Black" + (not boo) * "Blue"



def x2px(xval):
    """Convert an x value to a pixel value"""
    return zero_loc[0]+xval*px


def y2px(yval):
    """Convert an x value to a pixel value"""
    return zero_loc[1]-yval*px


def xytopx(coord):
    """If you give me a coordinate, I'll convert it
    from the real system to pixels."""
    return(zero_loc[0]+coord[0]*px, zero_loc[1]-coord[1]*px)


def draw_XB(canvas, XB, head_xy):
    """Draw the XB on the canvas for the first time
    """
    if canvas.find_withtag('fil') is not (): #if we've called before
        try:
            update_XB(canvas, XB, head_xy)       #knock on the other door
            return
        except TclError: # This is meant to keep the shell from hanging when
            return       # the window is force quit, but doesn't do so.
    #get xb locations
    [c_loc, h_loc] = XB.find_locs()
    #scale and translate locs
    c_loc = xytopx(c_loc)
    h_loc = xytopx(h_loc)
    head = xytopx(head_xy)
    #draw lines
    fil_id = canvas.create_line(x2px(-3), y2px(0),
                                x2px(15), y2px(0), width=4,
                                tag ="fil")
    neck_id = canvas.create_line(x2px(0), y2px(0),
                                 c_loc[0],c_loc[1],
                                 fill=col(XB.SPR_N), width=sw,
                                 tag="n")
    glob_id = canvas.create_line(c_loc[0],c_loc[1],
                                 h_loc[0],h_loc[1],
                                 fill=col(XB.SPR_G), width=sw,
                                 tag="g")
    thic_id = canvas.create_oval(x2px(0)-sr, y2px(0)-sr,
                                 x2px(0)+sr, y2px(0)+sr,
                                 fill=col(XB.SPR_T), outline=col(XB.SPR_T),
                                 tag="t")
    conv_id = canvas.create_oval(c_loc[0]-sr, c_loc[1]-sr,
                                 c_loc[0]+sr, c_loc[1]+sr,
                                 fill=col(XB.SPR_C), outline=col(XB.SPR_C),
                                 tag="c")
    head_id = canvas.create_rectangle(head[0]-sr, head[1]-sr,
                                      head[0]+sr, head[1]+sr,
                                      tag="head")
    #get and draw status text at top of canvas
    xb_text = 'XB Class: '+XB.give_me_curr_xb_class()
    bad_text = 'Badness: '+'%.3f' % XB.give_me_badness(head_xy)
    xb_id = canvas.create_text(10, 5, anchor=NW, text=xb_text, tag="xb")
    bad_id = canvas.create_text(10, 20, anchor=NW, text=bad_text, tag="bad")
    canvas.update()
    return (thic_id, neck_id, conv_id, glob_id, head_id)


def update_XB(canvas, XB, head_xy):
    #get xb locations
    [c_loc, h_loc] = XB.find_locs()
    #scale and translate locs
    c_loc = (zero_loc[0]+c_loc[0]*px, zero_loc[1]-c_loc[1]*px)
    h_loc = (zero_loc[0]+h_loc[0]*px, zero_loc[1]-h_loc[1]*px)
    head = xytopx(head_xy)
    #draw lines
    canvas.coords('n',zero_loc[0],zero_loc[1],c_loc[0],c_loc[1])
    canvas.coords('c',c_loc[0]-sr,c_loc[1]-sr,c_loc[0]+sr,c_loc[1]+sr)
    canvas.coords('g',c_loc[0],c_loc[1],h_loc[0],h_loc[1])
    canvas.coords('head',head[0]-sr,head[1]-sr,head[0]+sr,head[1]+sr)
    #get and draw status text at top of canvas
    xb_text = 'XB Class: '+XB.give_me_curr_xb_class()
    bad_text = 'Badness: '+'%.2f' % XB.give_me_badness(head_xy)
    canvas.itemconfigure("xb", text=xb_text)
    canvas.itemconfigure("bad", text=bad_text)
    canvas.update()


def new_root():
    #Make a new root window and name it
    root = Tk()
    root.title("XB")
    root.resizable(0, 0)
    #Create a frame to hold the canvas we will draw on
    frame = Frame(root, bd=3)
    frame.pack()
    #Create our canvas
    canvas = Canvas(frame, width=200, height=250, bd=0, highlightthickness=0)
    canvas.pack()
    #Display the whole thing
    root.update()
    return canvas


def kill_root(canvas):
    #destroy the root window containing a canvas of interest
    canvas.winfo_toplevel().destroy()


#canvas = new_root() #now done from calling environment

#use root.update_idletasks() to redraw root window
#use root.update() to process events such as window closings
