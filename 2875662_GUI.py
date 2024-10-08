#Importing necessary modules
from tkinter import *
from tkinter import Button
import tkinter as tk
from tkinter import ttk
import sqlite3
from PIL import Image, ImageTk
import PIL.Image
import io
from tkinter import filedialog
import csv

#Create the root window
root= Tk()
root.title("Chemical Compounds") #Set title
root.configure(bg= "snow2") #Set background colour
root.geometry("500x500") #Set window size

#Define the style for the treeview widget
style = ttk.Style()

style.configure("Treeview",
    background= "#D3D3D3",
    foreground= "black",
    rowheight= 100,
    fieldbackground= "#D3D3D3")

style.map("Treeview",
    background= [('selected', "#347083")])

style.configure("Treeview.Heading", font=("Arial", 10, "bold"))

#Create a frame for the table
table_frame = Frame(root)
table_frame.place(x=350, y=50, width=1080, height=760)  

#Create a frame for filters
data_frame = LabelFrame(root, text="Filters", width=300)
data_frame.config(font=("Arial", 12, "bold"))
data_frame.place(x=10, y=50, width=330, height=760) 

#Create a label for displaying the number of compounds
compounds_label = Label(data_frame, text="Number of Compounds: 0", font=("Arial", 10, "bold"))
compounds_label.grid(row=16, column=0, columnspan=2, padx=10, pady=10)

#Addong a vertical scrollbar
tree_scroll = Scrollbar(table_frame)
tree_scroll.pack(side=RIGHT, fill=Y)

#Adding a horizontal scrollbar
tree_scroll_x = Scrollbar(table_frame, orient=HORIZONTAL)
tree_scroll_x.pack(side=BOTTOM, fill=X)


# Create Treeview widget
tree = ttk.Treeview(table_frame, yscrollcommand=tree_scroll.set, xscrollcommand=tree_scroll_x.set, selectmode="extended")
tree.pack(expand=True, fill="both")

# Configure the scrollbars to control the Treeview widget
tree_scroll.config(command=tree.yview)
tree_scroll_x.config(command=tree.xview)

#Creating column names
tree["columns"] = ("ID", "Name", "Mass", "LogP", "LogD", "Rotatable Bonds", "Ring Count", "Donor Count", "Acceptor Count", "PSA", "FAR" , "Refractivity", "Total Atoms")

#format tree
tree.column("#0", anchor= CENTER, width= 120)
tree.column("ID", anchor= CENTER, width=80)
tree.column("Mass", anchor= CENTER, width= 80)
tree.column("Name", anchor= W, width= 100)
tree.column("LogP", anchor= CENTER, width= 80)
tree.column("LogD", anchor= CENTER, width= 80)
tree.column("Rotatable Bonds", anchor= CENTER, width= 100)
tree.column("Ring Count", anchor= CENTER, width= 100)
tree.column("Donor Count", anchor= CENTER, width= 100)
tree.column("Acceptor Count", anchor= CENTER, width= 110)
tree.column("PSA", anchor= CENTER, width= 80)
tree.column("FAR", anchor= CENTER, width= 80)
tree.column("Refractivity",anchor= CENTER, width= 100)
tree.column("Total Atoms", anchor= CENTER, width= 100)

tree.heading("#0", text= "Structure", anchor = CENTER)
tree.heading("ID", text="ID", anchor= CENTER )
tree.heading("Name", text="Name", anchor= W)
tree.heading("Mass", text= "Mass", anchor= CENTER)
tree.heading("LogP", text="LogP", anchor= CENTER)
tree.heading("LogD", text="LogD", anchor= CENTER)
tree.heading("Rotatable Bonds", text="Rotatable Bonds", anchor= CENTER)
tree.heading("Ring Count", text="Ring Count", anchor= CENTER)
tree.heading("Donor Count", text= "Donor Count" , anchor= CENTER)
tree.heading("Acceptor Count",text= "Acceptor Count", anchor= W)
tree.heading("PSA", text= "PSA", anchor= CENTER)
tree.heading("FAR", text= "FAR", anchor= CENTER)
tree.heading("Refractivity", text= "Refractivity", anchor = CENTER)
tree.heading("Total Atoms", text= "Total Atoms", anchor= CENTER)

tree.tag_configure('oddrow', background='azure2')
tree.tag_configure('evenrow', background='aliceblue')

# Insert fetched data into Treeview
image_list= []
# Define a function to populate the Treeview with data from the database
def populate_tree():
    # Connect to the SQLite database
    conn = sqlite3.connect("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/drugs_python/chem_db")
    c = conn.cursor()

    # Clear existing rows in the Treeview
    for row in tree.get_children():
        tree.delete(row)

    # Execute a SQL query to fetch all data from the molecules table
    c.execute("SELECT id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, FAR, refractivity, total_atoms, image FROM molecules")
    rows = c.fetchall()

    # Counter for alternating row colors
    count = 0
    # Loop through each row fetched from the database
    for row in rows:
        # Unpack the row data into variables
        id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms,img= row  
        # Open the image from binary data
        img = Image.open(io.BytesIO(img))
        img = img.resize((100,100), PIL.Image.BILINEAR)
        img = ImageTk.PhotoImage(img)
        image_list.append(img)
        
        # Insert the row data into the Treeview
        if count % 2 == 0:
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms), tags= ('evenrow',) )
        else :
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms), tags= ('oddrow',) )  
        count +=1  

    conn.close()
    compounds_label.config(text=f"Number of Compounds: {len(rows)}")

# Call the populate_tree function to initially populate the Treeview
populate_tree()

# Define a function to create a label and entry box in the data frame
def create_label_and_entry(parent, row, label_text, entry_width):
    label = Label(parent, text=label_text)
    label.grid(row=row, column=0, padx=10, pady=10)
    entry = Entry(parent, width=entry_width)
    entry.grid(row=row, column=1, padx=10, pady=10)
    return entry

# Create labels and entry boxes
mass_entry = create_label_and_entry(data_frame, 0, "Mass", 10)
logp_entry = create_label_and_entry(data_frame, 1, "LogP", 10)
logd_entry = create_label_and_entry(data_frame, 2, "LogD", 10)
donor_count_entry = create_label_and_entry(data_frame, 3, "Donor Count", 10)
acc_count_entry = create_label_and_entry(data_frame, 4, "Acceptor Count", 10)
ring_count_entry = create_label_and_entry(data_frame, 5, "Ring Count", 10)
FAR_entry = create_label_and_entry(data_frame, 6, "FAR", 10)
PSA_entry = create_label_and_entry(data_frame, 7, "PSA", 10)
refrac_entry = create_label_and_entry(data_frame, 8, "Refractivity", 10)
atoms_entry= create_label_and_entry(data_frame, 9, "Total Atoms", 10)
rot_entry= create_label_and_entry(data_frame, 10, "Rotatable Bonds", 10)

# Define a function to add an entry box with a "to" label

def add_entry_with_to_label(row, text, width):
    label = Label(data_frame, text=text)
    label.grid(row=row, column=2, padx=5, pady=10)
    entry = Entry(data_frame, width= width)
    entry.grid(row=row, column=3, padx=5, pady=10)
    return entry

# Add entries with "to" labels
to_mass = add_entry_with_to_label(0, "to", 10)
to_logp = add_entry_with_to_label(1, "to", 10)
to_logd = add_entry_with_to_label(2, "to", 10)
to_donor_count = add_entry_with_to_label(3, "to", 10)
to_acc_count = add_entry_with_to_label(4, "to", 10)
to_ring_count = add_entry_with_to_label(5, "to", 10)
to_FAR = add_entry_with_to_label(6, "to", 10)
to_PSA = add_entry_with_to_label(7, "to", 10)
to_refrac= add_entry_with_to_label(8, "to", 10)
to_atoms= add_entry_with_to_label(9, "to", 10)
to_rot= add_entry_with_to_label(10, "to", 10)

# Define a function to apply custom filters based on user input
def apply_filter():
    # Establish connection to the database
    conn = sqlite3.connect("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/drugs_python/chem_db")
    c = conn.cursor()

    # Define the base query
    query = '''SELECT id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, FAR, refractivity, total_atoms, image
               FROM molecules 
               WHERE 1=1'''

    # Define parameters and their corresponding entry fields
    params = {
        "mol_weight": (mass_entry.get(), to_mass.get()),
        "logp": (logp_entry.get(), to_logp.get()),
        "logd": (logd_entry.get(), to_logd.get()),
        "donor_count": (donor_count_entry.get(), to_donor_count.get()),
        "acceptor_count": (acc_count_entry.get(), to_acc_count.get()),
        "ring_count": (ring_count_entry.get(), to_ring_count.get()),
        "FAR": (FAR_entry.get(), to_FAR.get()),
        "psa": (PSA_entry.get(), to_PSA.get()),
        "refractivity": (refrac_entry.get(), to_refrac.get()),
        "total_atoms": (atoms_entry.get(), to_atoms.get()),
        "rotatable_bonds": (rot_entry.get(), to_rot.get())
    }

    # Add conditions to the query based on provided input
    for param, (from_value, to_value) in params.items():
        if from_value and to_value:
            query += f" AND {param} BETWEEN {from_value} AND {to_value}"

    # Execute the query
    c.execute(query)

    # Fetch the filtered data
    rows = c.fetchall()

    #Updating the number of compounds based on filtered results
    compounds_label.config(text=f"Number of Compounds: {len(rows)}")
    # Clear existing entries in the Treeview
    tree.delete(*tree.get_children())
    image_list.clear()
    count = 0
    # Insert filtered data into the Treeview
    for row in rows:
        id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms, img= row  # Assuming the image is the last column in the table
        img = Image.open(io.BytesIO(img))
        img = img.resize((100,100), PIL.Image.BILINEAR)
        img = ImageTk.PhotoImage(img)
        image_list.append(img)
        
        if count % 2 == 0:
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms), tags= ('evenrow',) )
        else :
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms), tags= ('oddrow',) )  
        count +=1  

    conn.close()
    
#Button for applying custom filters
filter_button = Button(data_frame, text="Apply Custom Filters", command=apply_filter, bg= "white")
filter_button.grid(row=11, column=0, columnspan=2, padx=10, pady=10)

# Define a function to clear all parameter entry boxes
def clear_parameters():
    # Clear all entry boxes
    mass_entry.delete(0, END)
    to_mass.delete(0, END)
    logp_entry.delete(0, END)
    to_logp.delete(0, END)
    logd_entry.delete(0, END)
    to_logd.delete(0, END)
    donor_count_entry.delete(0, END)
    to_donor_count.delete(0, END)
    acc_count_entry.delete(0, END)
    to_acc_count.delete(0, END)
    ring_count_entry.delete(0, END)
    to_ring_count.delete(0, END)
    FAR_entry.delete(0, END)
    to_FAR.delete(0, END)
    PSA_entry.delete(0, END)
    to_PSA.delete(0, END)
    refrac_entry.delete(0, END)
    to_refrac.delete(0, END)
    atoms_entry.delete(0, END)
    to_atoms.delete(0, END)
    rot_entry.delete(0, END)
    to_rot.delete(0, END)

# Button to clear all parameters
clear_params_button = Button(data_frame, text="Clear Parameters", command=clear_parameters, bg= "white")
clear_params_button.grid(row=11, column=2, columnspan=2, padx=10, pady=10)


#Function to filter data based on lipinski's rule
def filter_lipinski_rule():
    
    # Fetch data from the database satisfying the Lipinski rule
    conn = sqlite3.connect("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/drugs_python/chem_db")
    c = conn.cursor()

    #Fetching rows from the database where lipinski's rule is true
    query= '''SELECT id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, FAR, refractivity, total_atoms, image FROM molecules WHERE lipinski_rule = 1'''
    c.execute(query)
    rows = c.fetchall()

    #Updating the number of compounds based on the filtered results
    compounds_label.config(text=f"Number of Compounds: {len(rows)}")

    #Setting predefined entries to be displayed for the cut off values when filter applied 
    mass_entry.delete(0, END)
    mass_entry.insert(0, "")
    to_mass.delete(0, END)
    to_mass.insert(0, "500")
    to_logp.delete(0, END)
    to_logp.insert(0, "5")
    donor_count_entry.delete(0, END)
    donor_count_entry.insert(0, " ")
    to_donor_count.delete(0, END)
    to_donor_count.insert(0, "5")
    acc_count_entry.delete(0, END)
    acc_count_entry.insert(0, " ")
    to_acc_count.delete(0, END)
    to_acc_count.insert(0, "10")

   

    # Clear existing entries in the Treeview
    tree.delete(*tree.get_children())
    image_list.clear()
    count = 0
    for row in rows:
        id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms, img= row  # Assuming the image is the last column in the table
        img = Image.open(io.BytesIO(img))
        img = img.resize((100,100), PIL.Image.BILINEAR)
        img = ImageTk.PhotoImage(img)
        image_list.append(img)
        
        if count % 2 == 0:
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count , refractivity, total_atoms), tags= ('evenrow',) )
        else :
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count , refractivity, total_atoms), tags= ('oddrow',) )  
        count +=1  

    conn.close()
     # Apply the Lipinski rule filter
    apply_filter()
    
# Button to trigger Lipinski rule filtering
lipinski_button = Button(data_frame, text="Lipinski Rule", command=filter_lipinski_rule, bg= "white")
lipinski_button.grid(row=12, column=0, columnspan=2, padx=10, pady=10)

#Creating function for lead-likeness
def filter_lead_likeness():
    # Fetch data from the database satisfying the Lipinski rule
    conn = sqlite3.connect("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/drugs_python/chem_db")
    c = conn.cursor()
    query= '''SELECT id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, FAR, refractivity, total_atoms, image FROM molecules WHERE lead_likeness = 1'''
    c.execute(query)
    rows = c.fetchall()

    compounds_label.config(text=f"Number of Compounds: {len(rows)}")
    mass_entry.delete(0, END)
    mass_entry.insert(0, " ")
    to_mass.delete(0, END)
    to_mass.insert(0, "450")
    logd_entry.delete(0, END)
    logd_entry.insert(0, "-4")
    to_logd.delete(0, END)
    to_logd.insert(0, "4")
    donor_count_entry.delete(0, END)
    donor_count_entry.insert(0, " ")
    to_donor_count.delete(0, END)
    to_donor_count.insert(0, "5")
    acc_count_entry.delete(0, END)
    acc_count_entry.insert(0, " ")
    to_acc_count.delete(0, END)
    to_acc_count.insert(0, "8")
    ring_count_entry.delete(0, END)
    ring_count_entry.insert(0, " ")
    to_ring_count.delete(0, END)
    to_ring_count.insert(0, "4")
    rot_entry.delete(0, END)
    rot_entry.insert(0, " ")
    to_rot.delete(0, END)
    to_rot.insert(0, "10")
    
    
    
    tree.delete(*tree.get_children())
    image_list.clear()
    count = 0
    for row in rows:
        id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms, img= row  # Assuming the image is the last column in the table
        img = Image.open(io.BytesIO(img))
        img = img.resize((100,100), PIL.Image.BILINEAR)
        img = ImageTk.PhotoImage(img)
        image_list.append(img)
        
        if count % 2 == 0:
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count , refractivity, total_atoms), tags= ('evenrow',) )
        else :
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count , refractivity, total_atoms), tags= ('oddrow',) )  
        count +=1  

    conn.close()
    # Apply the Lipinski rule filter
    apply_filter()

# Button to trigger Lead-likeness filtering
lead_button = Button(data_frame, text="Lead Likeness", command=filter_lead_likeness, bg= "white")
lead_button.grid(row=12, column=2, columnspan=2, padx=10, pady=10)

#Creating a function for bioavailability
def filter_bioavail():
    # Fetch data from the database satisfying the Lipinski rule
    conn = sqlite3.connect("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/drugs_python/chem_db")
    c = conn.cursor()
    query= '''SELECT id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, FAR, refractivity, total_atoms, image FROM molecules WHERE bioavailability = 1'''
    c.execute(query)
    rows = c.fetchall()
    compounds_label.config(text=f"Number of Compounds: {len(rows)}")
    mass_entry.delete(0, END)
    mass_entry.insert(0, " ")
    to_mass.delete(0, END)
    to_mass.insert(0, "500")
    logp_entry.delete(0, END)
    logp_entry.insert(0, " ")
    to_logp.delete(0, END)
    to_logp.insert(0, "5")
    donor_count_entry.delete(0, END)
    donor_count_entry.insert(0, " ")
    to_donor_count.delete(0, END)
    to_donor_count.insert(0, "5")
    acc_count_entry.delete(0, END)
    acc_count_entry.insert(0, " ")
    to_acc_count.delete(0, END)
    to_acc_count.insert(0, "10")
    FAR_entry.delete(0, END)
    FAR_entry.insert(0, " ")
    to_FAR.delete(0, END)
    to_FAR.insert(0, "5")
    rot_entry.delete(0, END)
    rot_entry.insert(0, " ")
    to_rot.delete(0, END)
    to_rot.insert(0, "10")
    PSA_entry.delete(0, END)
    PSA_entry.insert(0, " ")
    to_PSA.delete(0, END)
    to_PSA.insert(0, "200")
        
    

    # Clear existing entries in the Treeview
    tree.delete(*tree.get_children())
    image_list.clear()
    count = 0
    for row in rows:
        id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms, img= row  # Assuming the image is the last column in the table
        img = Image.open(io.BytesIO(img))
        img = img.resize((100,100), PIL.Image.BILINEAR)
        img = ImageTk.PhotoImage(img)
        image_list.append(img)
        
        if count % 2 == 0:
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count , refractivity, total_atoms), tags= ('evenrow',) )
        else :
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms), tags= ('oddrow',) )  
        count +=1  

    conn.close()
    # Apply the Lipinski rule filter
    apply_filter()
    

# Button to trigger bioavailability rule filtering
bioavail_button = Button(data_frame, text="Bioavailability", command=filter_bioavail, bg= "white")
bioavail_button.grid(row=13, column=0, columnspan=2, padx=10, pady=10)

#Creating a function for ghose
def filter_ghose():
    # Fetch data from the database satisfying the Lipinski rule
    conn = sqlite3.connect("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/drugs_python/chem_db")
    c = conn.cursor()
    query= '''SELECT id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, FAR, refractivity, total_atoms, image FROM molecules WHERE ghose = 1'''
    c.execute(query)
    rows = c.fetchall()
    compounds_label.config(text=f"Number of Compounds: {len(rows)}")
    mass_entry.delete(0, END)
    mass_entry.insert(0, "160")
    to_mass.delete(0, END)
    to_mass.insert(0, "480")
    logp_entry.delete(0, END)
    logp_entry.insert(0, "-0.4")
    to_logp.delete(0, END)
    to_logp.insert(0, "5.6")
    refrac_entry.delete(0, END)
    refrac_entry.insert(0, "40")
    to_refrac.delete(0, END)
    to_refrac.insert(0, "130")
    atoms_entry.delete(0, END)
    atoms_entry.insert(0, "20")
    to_atoms.delete(0, END)
    to_atoms.insert(0, "70")

    
    # Clear existing entries in the Treeview
    tree.delete(*tree.get_children())
    image_list.clear()
    count = 0
    for row in rows:
        id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms ,img= row  # Assuming the image is the last column in the table
        img = Image.open(io.BytesIO(img))
        img = img.resize((100,100), PIL.Image.BILINEAR)
        img = ImageTk.PhotoImage(img)
        image_list.append(img)
        
        if count % 2 == 0:
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms), tags= ('evenrow',) )
        else :
            tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms), tags= ('oddrow',) )  
        count +=1  

    conn.close()
    # Apply the Lipinski rule filter
    apply_filter()

    

# Button to trigger ghose rule filtering
ghose_button = Button(data_frame, text="Ghose", command=filter_ghose, bg= "white")
ghose_button.grid(row=13, column=2, columnspan=2, padx=10, pady=10)

reset_button = Button(data_frame, text="Reset", bg= "white")
reset_button.grid(row=14, column= 1, columnspan=2, padx=10, pady=10)
reset_button.config(command= populate_tree)

#Creating a function for search bar to searcg compound by name
def search():
    query = search_entry.get()
    if query:
        # Clear existing entries in the Treeview
        tree.delete(*tree.get_children())
        conn = sqlite3.connect("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/drugs_python/chem_db")
        c = conn.cursor()
        query1 = '''SELECT id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, FAR, refractivity, total_atoms, image FROM molecules WHERE name LIKE ?'''
        c.execute(query1, ('%' + query + '%',)) 
        rows = c.fetchall()
        compounds_label.config(text=f"Number of Compounds: {len(rows)}")

        image_list.clear()
        count = 0
        for row in rows:
            id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count, refractivity, total_atoms, img= row  # Assuming the image is the last column in the table
            img = Image.open(io.BytesIO(img))
            img = img.resize((120,120), PIL.Image.BILINEAR)
            img = ImageTk.PhotoImage(img)
            image_list.append(img)
            
            if count % 2 == 0:
                tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count , refractivity, total_atoms), tags= ('evenrow',) )
            else :
                tree.insert(parent= '', index= 'end', iid= count, image=img, values= (id, name, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, fused_aromatic_ring_count , refractivity, total_atoms), tags= ('oddrow',) )  
            count +=1 
       
        conn.close()
        
# Create a search bar
search_frame = Frame(root, bg="white")
search_frame.place(x=600, y=10, width=320, height=30)
search_entry = Entry(search_frame, font=('Arial', 10), bg="white", relief=GROOVE)
search_entry.pack(side=LEFT, expand=True, fill="both", padx=5, pady=5)
search_button = Button(search_frame, text="Search Compound", command=search)
search_button.pack(side=LEFT, padx=5, pady=5)

# Bind the <Return> event to the search function
search_entry.bind("<Return>", lambda event: search())

def save_filtered_data():
    # Get filtered data from the Treeview
    filtered_data = []
    for item in tree.get_children():
        values = tree.item(item)['values']
        filtered_data.append(values)
    
    # Ask user to select file path for saving
    file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
    
    # Write data to the file
    with open(file_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["ID", "Name", "Mass", "LogP", "LogD", "Rotatable Bonds", "Ring Count", "Donor Count", "Acceptor Count", "PSA", "FAR", "Refractivity", "Total Atoms"])
        for row in filtered_data:
            writer.writerow(row)

# Button to save filtered data
save_button = Button(data_frame, text="Save Filtered data", command=save_filtered_data, bg= "white")
save_button.grid(row=15, column=1, columnspan=2, padx=10, pady=10)

tree.pack(expand=True, fill="both")

root.mainloop()
