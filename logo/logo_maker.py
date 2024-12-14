import numpy as np
from scatter_letters import sl
from matplotlib import pyplot as plt
from numpy.linalg import norm


def detect_words_width(words, verbose = False):
    words_is_list = True
    if np.shape(words) == (): # single word
        words_is_list = False
        words = [words]

    width_words = []
    for cpt_word, myword in enumerate(words):
        mytext = sl.text_to_data(myword, intensity = 50, repeat=False)
        total_width_word = 0
        for cpt_letter, letter in enumerate(myword[:-1]):
            data_letter = np.array(mytext[cpt_letter])
            data_letter[0] -= np.min(data_letter[0])
            width_letter = np.max(data_letter[0])
            total_width_word += width_letter
        width_words.append(total_width_word)
    if verbose:
        print("Words width detected")
    if words_is_list:
        return width_words
    else:
        return width_words[0]
    
def place_the_letters(words, width_words, space_letter_x = 100,
                      space_letter_y = 100, space_box_y = 200,
                      verbose = False, show = True):
    n_letter = 0
    for cpt_word, myword in enumerate(words):
        mytext = sl.text_to_data(myword, intensity = 250, repeat=False)
        width_letter_total = 0
        for cpt_letter, letter in enumerate(mytext[:-1]):
            data_letter = np.array(mytext[cpt_letter])
            data_letter[0] -= np.min(data_letter[0])
            width_letter = np.max(data_letter[0])
            height_letter = np.max(data_letter[1])
            data_letter[0] += width_letter_total
            data_letter[0] += space_letter_x*cpt_letter
            data_letter[0] -= width_words[cpt_word]
            data_letter[0] -= (len(myword)-2)*space_letter_x
            data_letter[1] -= np.min(data_letter[1])
            data_letter[1] -= space_letter_y*cpt_word
            data_letter[1] -= height_letter*cpt_word
            width_letter_total += width_letter  
            if n_letter == 0:
                all_letters = data_letter
            else:
                all_letters = np.vstack([all_letters.T, data_letter.T]).T
            n_letter += 1
    all_letters = all_letters.astype('float64')
    all_letters[0] -= (np.max(all_letters[0]) + np.min(all_letters[0]))/2
    all_letters[1] -= (np.max(all_letters[1]) + np.min(all_letters[1]))/2       
    Lx = np.max(all_letters[0]) - np.min(all_letters[0]) + 2*space_letter_x
    Ly = np.max(all_letters[1]) - np.min(all_letters[1]) + 2*space_box_y
    if show:
        fig, ax = plt.subplots(figsize=(14,5))
        ax.plot(all_letters[0], all_letters[1], '.', color="b")
        ax.axis('equal')
        fig.tight_layout()
        plt.axis('off')
        plt.show()
    if verbose:
        print("Letters placed in space")
    return all_letters, Lx, Ly

def place_the_atoms(all_letters, Lx, Ly, random_placement = True,
                    number_of_try = 3500, d0 = 10,
                    show = False, rescaling = 4.5):
    box = np.array([Lx, Ly])
    atoms = []
    cpt_atoms = 0
    if random_placement:
        for N in range(number_of_try):
            cpt_atoms += 1
            x = np.random.random()*Lx - Lx/2
            y = np.random.random()*Ly - Ly/2
            d = np.min(norm(np.remainder(all_letters.T - np.array([x,y]) + box/2., box) - box/2., axis=1))
            if d < 1:
                atoms.append([cpt_atoms, 1, x/rescaling, y/rescaling, 0])
            else:
                atoms.append([cpt_atoms, 2, x/rescaling, y/rescaling, 0])
    else:
        shift_x = 0
        shift_y = 0

        if shift_x % 2 == 0:
            x = - Lx/2 +d0/2
        else:
            x = - Lx/2    

        while x < Lx/2:
            if shift_y % 2 == 0:
                y = -Ly/2+d0/2
            else:
                y = -Ly/2

            while y < Ly/2:
                cpt_atoms += 1
                d = np.min(norm(np.remainder(all_letters.T - np.array([x,y]) + box/2., box) - box/2., axis=1))
                if d < 12:
                    atoms.append([cpt_atoms, 1, x/rescaling, y/rescaling, 0])
                else:
                    atoms.append([cpt_atoms, 2, x/rescaling, y/rescaling, 0])
                y += d0
            x += d0/2
            shift_y += 1

    
    atoms =  np.array(atoms)
    if show:
        fig, ax = plt.subplots(figsize=(14,5))
        ax.plot(atoms[atoms.T[1] == 2].T[2],
                atoms[atoms.T[1] == 2].T[3], 'o', color="b")
        ax.plot(atoms[atoms.T[1] == 1].T[2],
                atoms[atoms.T[1] == 1].T[3], 'o', color="y")
        ax.axis('equal')
        fig.tight_layout()
        plt.axis('off')
        plt.show()
    Lx /= rescaling
    Ly /= rescaling
    return atoms, Lx, Ly

def write_lammps_data(filename, atoms, Lx, Ly, increased_size = 5):
    f = open(filename, "w")
    f.write('# LAMMPS data file \n\n')
    f.write(str(len(atoms))+' atoms\n')
    f.write('\n')
    f.write('2 atom types\n')
    f.write('\n')
    f.write(str(-Lx/2-increased_size/2)+' '+str(Lx/2+increased_size/2)+' xlo xhi\n')
    f.write(str(-Ly/2-increased_size/2)+' '+str(Ly/2+increased_size/2)+' ylo yhi\n')
    f.write('-1 1 zlo zhi\n')
    f.write('\n')
    f.write('Atoms\n')
    f.write('\n')
    for nlin in range(len(atoms)):
        newline = atoms[nlin]
        for col in range(len(newline)):
            if col < 3:
                f.write(str(int(newline[col]))+' ')
            else :
                f.write(str(newline[col])+' ')
        f.write('\n')
    f.close()