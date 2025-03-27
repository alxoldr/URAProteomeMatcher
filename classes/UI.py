#!/usr/bin/python3.13

import json

from tkinter import *
from tkinter.filedialog import askopenfilename

def create_label(frame, text):
    """ Creates a label """
    label = Label(frame, text=text)
    label.pack(fill=X)
    return label

def create_entry(frame, text):
    """ Creates a text entry """
    entry = Entry(frame)
    entry.insert(0, text)
    entry.pack(fill=X)
    return entry

def create_button(frame, text, function, vars):
    """ Creates a button that can trigger a specified function """
    button = Button(frame, text=text, command=lambda:function(vars))
    button.pack(fill=X)
    return button

def create_optionmenu(frame, options):
    """ Creates an option menu """
    variable = StringVar(value=options[0])
    optionmenu = OptionMenu(frame, variable, *options)
    optionmenu.pack(anchor=W, padx=5, pady=5)
    return variable

def create_checkbox(frame, text):
    """ Creates a textbox """
    variable = BooleanVar(value=False)
    checkbox = Checkbutton(frame, text=text, variable=variable)
    checkbox.pack(anchor=W, padx=5, pady=5)
    return variable

def open_popup(root, title, text):
    """ Opens a popup window with a given message """
    popup = Toplevel(root)
    popup.title(title)

    label = Label(popup, text=text)
    label.pack(pady=20)

def get_file_path(entry):
    """ Prompts for a user to select a file, returns the file path  """
    Tk().withdraw()
    file_path = askopenfilename()
    entry.delete(0, END)
    entry.insert(0, file_path)
    return file_path
