import sqlite3
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Draw
import io

conn = sqlite3.connect("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/drugs_python/chem_db")
c = conn.cursor()

c.execute('''CREATE TABLE IF NOT EXISTS molecules
             (id INTEGER PRIMARY KEY,
             name TEXT,
             smiles TEXT,
             formula TEXT,
             mol_weight REAL,
             logp REAL,
             logd REAL,
             rotatable_bonds INTEGER,
             FAR INTEGER,
             ring_count INTEGER,
             donor_count INTEGER,
             acceptor_count INTEGER,
             psa REAL,
             refractivity INTEGER,
             total_atoms INTEGER,
             lipinski_rule BOOLEAN,
             bioavailability BOOLEAN,
             lead_likeness BOOLEAN,
             ghose BOOLEAN,
             image BLOB)''')

def calculate_lipinski_rule(mol):
    num_h_donors = Descriptors.NumHDonors(mol)
    num_h_acceptors = Descriptors.NumHAcceptors(mol)
    mol_weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)

    if num_h_donors <= 5 and num_h_acceptors <= 10 and mol_weight <= 500 and logp <= 5:
        return True
    else:
        return False
    
def calculate_lead_likeness(mol):
    num_h_donors = Descriptors.NumHDonors(mol)
    num_h_acceptors = Descriptors.NumHAcceptors(mol)
    mol_weight = Descriptors.MolWt(mol)
    logd = float(mol.GetProp('LogD'))
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    ring_count = Descriptors.RingCount(mol)
    if num_h_donors <= 5 and num_h_acceptors <= 8 and mol_weight <= 450 and -4 <= logd <= 4 and ring_count <= 4 and rotatable_bonds <=10:
        return True
    else:
        return False
    
def calculate_bioavailability(mol, fused_aromatic_ring):
    mol_weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    donor_count = Descriptors.NumHDonors(mol)
    acceptor_count = Descriptors.NumHAcceptors(mol)
    psa = Descriptors.TPSA(mol)  # Calculate PSA
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)

    if psa <= 200 and mol_weight <= 500 and logp <= 5 and donor_count <= 5 and acceptor_count <= 10 and rotatable_bonds <=10 and fused_aromatic_ring <=5:
        return True
    else:
        return False
    
def ghose(mol):
    mol_weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    refractivity = Descriptors.MolMR(mol)
    total_atoms= mol.GetNumAtoms()

    if -0.4 <= logp <= 5.6 and 160 <= mol_weight <= 480 and 40 <= refractivity <= 130 and 20 <= total_atoms <= 70:
        return True
    else:
        return False
        

molecules = []
for idx, mol in enumerate(Chem.SDMolSupplier("C:/Users/Aakanksha Choudhary/OneDrive/Desktop/drugs_python/Molecules4.sdf"), start=1):
    print(idx)
    if mol:
        name = mol.GetProp('Name') if mol.HasProp('Name') else ''
        formula = mol.GetProp('Formula') if mol.HasProp('Formula') else ''
        mol_weight = float(mol.GetProp('MolWeight')) if mol.HasProp('MolWeight') else 0.0
        smiles = Chem.MolToSmiles(mol)
        logp = round(Descriptors.MolLogP(mol), 3)
        logd = round(float(mol.GetProp('LogD')), 3)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        ring_count = Descriptors.RingCount(mol)
        donor_count = Descriptors.NumHDonors(mol)
        acceptor_count = Descriptors.NumHAcceptors(mol)
        psa = round(Descriptors.TPSA(mol), 3)
        refractivity = round(Descriptors.MolMR(mol), 3)
        total_atoms= mol.GetNumAtoms()

        rotatable_bonds_indexes = []

                
        rings = mol.GetRingInfo()
        aro_rings = mol.GetRingInfo().AtomRings()
        fused_aromatic_ring_count = 0
        aromatic_rings = []

        for ring_idx in range(rings.NumRings()):
            # Check if the ring is fused
            if rings.IsRingFused(ring_idx):
                # Check if the ring is aromatic
                is_aromatic = False
                for atom_idx in rings.AtomRings()[ring_idx]:
                    if mol.GetAtomWithIdx(atom_idx).GetIsAromatic():
                        is_aromatic = True
                        break
                # Increment the count if it's a fused aromatic ring
                if is_aromatic:
                    fused_aromatic_ring_count += 1

        lipinski_rule = True if calculate_lipinski_rule(mol) else False
        bioavailability = True if calculate_bioavailability(mol, fused_aromatic_ring_count) else False
        lead_likeness = True if calculate_lead_likeness(mol) else False
        ghose_filter = True if ghose(mol) else False

        img = Draw.MolToImage(mol, size=(300, 300))
        img_byte_array = io.BytesIO()
        img.save(img_byte_array, format='PNG')
        img_byte_array = img_byte_array.getvalue()

        molecules.append((idx, name, smiles, formula, mol_weight, logp, logd, rotatable_bonds, ring_count, donor_count, acceptor_count, psa, img_byte_array))

        c.execute('''INSERT INTO molecules (id, name, smiles, formula, mol_weight, logp, logd, rotatable_bonds, FAR, ring_count, donor_count, acceptor_count, psa, refractivity, total_atoms, lipinski_rule, bioavailability, lead_likeness, ghose, image)
                     VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                  (idx, name, smiles, formula, mol_weight, logp, logd, rotatable_bonds, fused_aromatic_ring_count, ring_count, donor_count, acceptor_count, psa, refractivity, total_atoms, lipinski_rule, bioavailability, lead_likeness, ghose_filter, sqlite3.Binary(img_byte_array)))
        
conn.commit()
conn.close()