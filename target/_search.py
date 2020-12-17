"""
Input Checker
===============

Corrects user input such as 'wasp43b' to 'WASP-43 b'
"""

def _input_checker(target):
    # Check WASP entry
    if target.lower().__contains__('wasp'):
        WASP = target[0:4].upper() + '-'
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{WASP}{numbers[0]}{last_char}'
    # Check TOI entry
    if target.lower().__contains__('toi'):
        TOI = target[0:3].upper() + '-'
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{TOI}{numbers[0]}{last_char}'
    # Check LHS entry
    if target.lower().__contains__('lhs'):
        LHS = target[0:3].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{LHS}{numbers[0]}{last_char}'
    # Check LTT entry
    if target.lower().__contains__('ltt'):
        LTT = target[0:3].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{LTT}{numbers[0]}{last_char}'
    # Check GJ entry
    if target.lower().__contains__('gj'):
        GJ = target[0:2].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{GJ}{numbers[0]}{last_char}'
    # Check HD entry
    if target.lower().__contains__('hd'):
        HD = target[0:2].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{HD}{numbers[0]}{last_char}'
    # Check HR entry
    if target.lower().__contains__('hr'):
        HR = target[0:2].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{HR}{numbers[0]}{last_char}'
    # Check NGTS entry
    if target.lower().__contains__('ngts'):
        NGTS = target[0:4].upper() + '-'
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{NGTS}{numbers[0]}{last_char}'
    # Check LP entry
    if target.lower().__contains__('lp'):
        LP = target[0:2].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify)
        last_char = f' {target[-1]}'
        target =  f'{LP}{numbers[0:3]}-{numbers[3:5]}{last_char}'
    # Check WD entry
    if target.lower().__contains__('wd'):
        WD = target[0:2].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify)
        last_char = f' {target[-1]}'
        target =  f'{WD}{numbers[0:4]}+{numbers[4:7]}{last_char}'
    # Check L entry
    # if target.lower().__contains__('l'):
    #     L = target[0:1].upper() + ' '
    #     classify = list(map(lambda sub:str(''.join( 
    #                     [i for i in sub if i.isnumeric()])), target))
    #     numbers = ''.join(classify)
    #     last_char = f' {target[-1]}'
    #     target =  f'{L}{numbers[0:3]}-{numbers[3:4]}{last_char}'