#!/usr/bin/env python3
import os 

def get_potcar_data():
    home = os.path.expanduser('~')
    data_potcars_file = home + '/bin/q-robot/books/potpaw_PBE.54/data_potcars'
    file_in = open(data_potcars_file, 'r')
    data = file_in.read()
    file_in.close()
    return eval(data)    
# elements = get_potcar_data().keys()

def concatenate(ele_list):
    ''' cat POTCAR1 POTCAR2 >> POTCAR '''
    home = os.path.expanduser('~')
    potcar_path = home + '/bin/q-robot/books/potpaw_PBE.54'
    f_out = open('POTCAR', 'w')
    for i in ele_list:
        print('\nAdd %s to the POTCAR' %(i))
        try: 
            ele_path = potcar_path + '/' + i
            f = open(ele_path+'/POTCAR', 'r')
            lines = f.readlines()
            f.close()
        except: 
            ele_path = potcar_path + '/' + i + '_sv'
            f = open(ele_path+'/POTCAR', 'r')
            lines = f.readlines()
            f.close()
        f_out.writelines(lines)
    f_out.close()

def get_potcar_infor(potcar, num=7):    
    ''' Get the information for single POTCAR file '''
    dict_ele_potcar = {}
    f = open(potcar, 'r')
    lines = f.readlines()
    f.close()
    infor = lines[num-7].strip().split()
    dict_ele_potcar['TYPE'] = infor[0]
    dict_ele_potcar['NAME'] = infor[1]
    dict_ele_potcar['DATE'] = infor[2]
    
    state = lines[num-4].strip().split('=')
    dict_ele_potcar[state[0].strip()] = state[1]
    
    eatom = lines[num-2].strip().split('=')
    dict_ele_potcar[eatom[0].strip()] = eatom[1].split()[0]
    
    for line in lines[num+1:num+18]:
        if '=' in line: 
            if ';' in line:
                line_ele = line.rstrip().split(';')
                for i in line_ele:
                    item_ele = [y.strip() for y in i.split('=')]
                    dict_ele_potcar[item_ele[0]] = item_ele[1].split()[0]
            else:
                item_ele = line.rstrip().split('=')
                dict_ele_potcar[item_ele[0].strip()] =  item_ele[1].split()[0].strip()
    return dict_ele_potcar

def get_multiple_potcar_infor(potcar):
    '''Get the information for concatenated  POTCAR with several elements''' 
    dict_potcars = {}
    f = open(potcar, 'r')
    lines = f.readlines()
    f.close()
    ele_list = []
    for num, line in enumerate(lines):
        if 'TITEL' in line:
            pot_name = line.rstrip().split()[3]
            ele_list.append(pot_name)
            dict_ele_potcar  = get_potcar_infor(potcar, num)
            dict_potcars.update({pot_name:dict_ele_potcar})
    return dict_potcars, ele_list

def get_potcars_infor(version):    
    '''Generate the data_potcars file for Q-robot to remember.''' 
    'version is the suffix in the potcar folders  of  potpaw_PBE.52, so it should be 52 or 54'
    home = os.path.expanduser('~')
    potcar_path = home + '/bin/q-robot/books/potpaw_PBE.' + version + '/'
    elements = [f for f in os.listdir(potcar_path) if os.path.isdir(potcar_path + f)]
    dict_potcars = {}
    for i in elements:
        potcar = potcar_path + '/' + i + '/POTCAR'
        dict_i = get_potcar_infor(potcar)
        dict_potcars[i] = dict_i    
    dict_name = home + '/bin/q-robot/books/potpaw_PBE.' + version + '/data_potcars'
    file_out = open(dict_name, 'w')
    file_out.write(str(dict_potcars))
    return dict_potcars
### Only run it when you update your potcars. 
#get_potcars_infor('52')
#get_potcars_infor('54') 
    
def read_potcar(potcar):
    ''' Print out the basic information of POTCAR in the current folder'''
    dict_potcars, ele_list = get_multiple_potcar_infor(potcar)    
    print('\n' *2)
    print('%s POTCAR Information %s\n' %('%'*25, '%'*25) )
    print_list = ['NAME', 'DATE', 'ZVAL', 'ENMAX', 'POMASS', 'TYPE', 'VRHFIN']
    print('NAME\tDATE\t\tN_ELE\tENMAX\tMASS\tTYPE\tE_STATE')
    for i in ele_list:
        v = dict_potcars.get(i)
        print_out = [v.get(items) for items in print_list]
        print('\t'.join(print_out))  
    print('\n%s Good Luck! %s\n' %('%'*29, '%'*29) )    
#read_potcar('POTCAR')    
    
