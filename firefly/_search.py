"""
Input Checker
===============

Corrects user input such as 'wasp43b' to 'WASP-43 b'
"""

def _is_anagram(str1, str2):
    list_str1 = list(str1)
    list_str1.sort()
    list_str2 = list(str2)
    list_str2.sort()
    return (list_str1 == list_str2)


def _input_checker(target):
    # Check WASP entry
    # if target.lower().__contains__('wasp'):
    if _is_anagram(target[0:4].lower(), 'wasp'):
        target = target.replace(target[0:4],'wasp')
        WASP = target[0:4].upper() + '-'
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{WASP}{numbers[0]}{last_char}'
    # Check TOI entry
    if _is_anagram(target[0:3].lower(), 'toi'):
        target = target.replace(target[0:3],'toi')
        TOI = target[0:3].upper() + '-'
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{TOI}{numbers[0]}{last_char}'
    # Check LHS entry
    if _is_anagram(target[0:3].lower(), 'lhs'):
        target = target.replace(target[0:3],'lhs')
        LHS = target[0:3].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{LHS}{numbers[0]}{last_char}'
    # Check LTT entry
    if _is_anagram(target[0:3].lower(), 'ltt'):
        target = target.replace(target[0:3],'ltt')
        LTT = target[0:3].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{LTT}{numbers[0]}{last_char}'
    # Check GJ entry
    if _is_anagram(target[0:2].lower(), 'gj'):
        target = target.replace(target[0:2],'gj')
        GJ = target[0:2].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{GJ}{numbers[0]}{last_char}'
    # Check HD entry
    if _is_anagram(target[0:2].lower(), 'hd'):
        target = target.replace(target[0:2],'hd')
        HD = target[0:2].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{HD}{numbers[0]}{last_char}'
    # Check HR entry
    if _is_anagram(target[0:2].lower(), 'hr'):
        target = target.replace(target[0:2],'hr')
        HR = target[0:2].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{HR}{numbers[0]}{last_char}'
    # Check NGTS entry
    if _is_anagram(target[0:4].lower(), 'ngts'):
        target = target.replace(target[0:4],'ngts')
        NGTS = target[0:4].upper() + '-'
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify).split()
        last_char = f' {target[-1]}'
        return f'{NGTS}{numbers[0]}{last_char}'
    # Check LP entry
    if _is_anagram(target[0:2].lower(), 'lp'):
        target = target.replace(target[0:2],'lp')
        LP = target[0:2].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify)
        last_char = f' {target[-1]}'
        return  f'{LP}{numbers[0:3]}-{numbers[3:5]}{last_char}'
    # Check WD
    if _is_anagram(target[0:2].lower(), 'wd'):
        target = target.replace(target[0:2],'wd')
        WD = target[0:2].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify)
        last_char = f' {target[-1]}'
        return f'{WD}{numbers[0:4]}+{numbers[4:7]}{last_char}'
    # Check L entry
    if _is_anagram(target[0:1].lower(), 'l'):
        target = target.replace(target[0:1],'l')
        L = target[0:1].upper() + ' '
        classify = list(map(lambda sub:str(''.join( 
                        [i for i in sub if i.isnumeric()])), target))
        numbers = ''.join(classify)
        last_char = f' {target[-1]}'
        return  f'{L}{numbers[0:3]}-{numbers[3:4]}{last_char}'
    
