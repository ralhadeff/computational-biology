sed 's/ALA/A/ ; s/CYS/C/ ; s/ASP/D/ ; s/GLU/E/ ; s/PHE/F/ ; s/GLY/G/ ; s/HIS/H/ ; s/ILE/I/ ; s/LYS/K/ ; s/LEU/L/ ; s/MET/M/ ; s/ASN/N/ ; s/PRO/P/ ; s/GLN/Q/ ; s/ARG/R/; s/SER/S/ ; s/THR/T/ ; s/VAL/V/ ; s/TRP/W/ ; s/TYR/Y/ ; s/HIE/H/ ; s/CSS/C/ ; s/CYI/C/ ; s/.../X/' 

# this bash script is to be piped to convert grep or similar outputs from 3 letter code to 1 letter code

# full command for pdb to fasta: grep CA input.pdb | awk '{print $4}' | $c/amino_3_to_1.bash | tr '\n' 'Z' | sed 's/Z//g ; echo ""'
