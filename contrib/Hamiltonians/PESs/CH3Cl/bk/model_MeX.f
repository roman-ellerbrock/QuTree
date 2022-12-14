C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE FUNCS(N,Q,R,P,V,DV,NPAR,MA,PST,SKIP,Y)
C% 21 Septembre 2010 A Viel
C% changed bei T. Weike
C%
C% pst(*,  1): 'VEXEN:      ' => Exactation energy 1->e,2->a
C% pst(*,  2): 'TMC_CX:     ' => TMC for C-X distance
C% pst(*,  3): 'TMC_CH:     ' => TMC for C-H distance
C% pst(*,  4): 'BLANCK1:    ' 
C% pst(*,  5): 'BLANCK2:    '
C% pst(*,  6): 'BLANCK3:    '
C% pst(*,  7): 'BLANCK4:    '
C% pst(*,  8): 'BLANCK5:    '
C% pst(*,  9): 'BLANCK6:    '
C% pst(*, 10): 'BLANCK7:    '
C% pst(*, 11): '1AMODE:     ' =>  1| 1| 1| 1| 1| 1| 1| 1
C% pst(*, 12): '1BMODE:     ' =>  1| 1| 1| 1| 1| 1| 1| 1
C% pst(*, 13): '1CMODE:     ' =>  1| 1| 1| 1| 1| 1| 1| 1
C% pst(*, 14): '1VJT6EA:    ' =>  0| 1| 1| 1| 1| 2|  |
C% pst(*, 15): '1VJT6EB:    ' =>  0| 1| 1| 1| 1| 2|  |
C% pst(*, 16): '1VJT6EC:    ' =>  0| 1| 1| 1| 1| 2|  |  
C% pst(*, 17): '1VEEAB:     ' =>  0| 1| 2| 4| 8|13|  |  
C% pst(*, 18): '1VEEAC:     ' =>  0| 1| 2| 4| 8|13|  |  
C% pst(*, 19): '1VEEBC:     ' =>  0| 1| 2| 4| 8|13|  |  
C% pst(*, 20): '1VEAA:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 21): '1VEBA:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 22): '1VECA:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 23): '1VEAB:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 24): '1VEBB:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 25): '1VECB:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 26): '1VEAC:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 27): '1VEBC:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 28): '1VECC:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 29): '1VAB:       ' =>  0| 1| 2| 3| 4| 5|  |  
C% pst(*, 30): '1VAC:       ' =>  0| 1| 2| 3| 4| 5|  |  
C% pst(*, 31): '1VBC:       ' =>  0| 1| 2| 3| 4| 5|  |
C% pst(*, 32): '1VAAA:      ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*, 33): '1VEEE:      ' =>  0| 0| 1| 6|  |  |  |   
C% pst(*, 34): '1VAEAEB:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*, 35): '1VBEAEB:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*, 36): '1VCEAEB:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*, 37): '1VAEAEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*, 38): '1VBEAEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*, 39): '1VCEAEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*, 40): '1VAEBEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*, 41): '1VBEBEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*, 42): '1VCEBEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*, 43): '1VABEA:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*, 44): '1VACEA:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*, 45): '1VBCEA:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*, 46): '1VABEB:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*, 47): '1VACEB:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*, 48): '1VBCEB:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*, 49): '1VABEC:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*, 50): '1VACEC:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*, 51): '1VBCEC:     ' =>  0| 0| 0| 1|  |  |  |   
C% pst(*, 52): '1VBLANCK1:  ' =>   |  |  |  |  |  |  |   
C% pst(*, 53): '1VBLANCK2:  ' =>   |  |  |  |  |  |  |   
C% pst(*, 54): '1VBLANCK3:  ' =>   |  |  |  |  |  |  |   
C% pst(*, 55): '1VBLANCK4:  ' =>   |  |  |  |  |  |  |   
C% pst(*, 56): '1VBLANCK5:  ' =>   |  |  |  |  |  |  |   
C% pst(*, 57): '1VBLANCK6:  ' =>   |  |  |  |  |  |  |   
C% pst(*, 58): '1VBLANCK7:  ' =>   |  |  |  |  |  |  |   
C% pst(*, 59): '1VBLANCK8:  ' =>   |  |  |  |  |  |  |   
C% pst(*, 60): '1VBLANCK9:  ' =>   |  |  |  |  |  |  |   
C% pst(*, 61): '1VBLANCK10: ' =>   |  |  |  |  |  |  |   
C% pst(*, 62): '1VBLANCK11: ' =>   |  |  |  |  |  |  |   
C% pst(*, 63): '1VBLANCK12: ' =>   |  |  |  |  |  |  |   
C% pst(*, 64): '1VBLANCK13: ' =>   |  |  |  |  |  |  |   
C% pst(*, 65): '1VBLANCK14: ' =>   |  |  |  |  |  |  |   
C% pst(*, 66): '1VBLANCK15: ' =>   |  |  |  |  |  |  |   
C% pst(*, 67): '1VBLANCK16: ' =>   |  |  |  |  |  |  |   
C% pst(*, 68): '1VBLANCK17: ' =>   |  |  |  |  |  |  |   
C% pst(*, 69): '1VBLANCK18: ' =>   |  |  |  |  |  |  |   
C% pst(*, 70): '1VBLANCK19: ' =>   |  |  |  |  |  |  |   
C% pst(*, 71): '1VBLANCK20: ' =>   |  |  |  |  |  |  |   
C% pst(*, 72): '1VBLANCK21: ' =>   |  |  |  |  |  |  |   
C% pst(*, 73): '1VBLANCK22: ' =>   |  |  |  |  |  |  |   
C% pst(*, 74): '1VBLANCK23: ' =>   |  |  |  |  |  |  |   
C% pst(*, 75): '1VBLANCK24: ' =>   |  |  |  |  |  |  |   
C% pst(*, 76): '1VBLANCK25: ' =>   |  |  |  |  |  |  |   
C% pst(*, 77): '1VBLANCK26: ' =>   |  |  |  |  |  |  |   
C% pst(*, 78): '1VBLANCK27: ' =>   |  |  |  |  |  |  |   
C% pst(*, 79): '1VBLANCK28: ' =>   |  |  |  |  |  |  |   
C% pst(*, 80): '1VBLANCK29: ' =>   |  |  |  |  |  |  |   
C% pst(*, 81): '1VBLANCK30: ' =>   |  |  |  |  |  |  |  
C% pst(*, 82): '2AMODE:     ' =>  1| 1| 1| 1| 1| 1| 1| 1
C% pst(*, 83): '2BMODE:     ' =>  1| 1| 1| 1| 1| 1| 1| 1
C% pst(*, 84): '2CMODE:     ' =>  1| 1| 1| 1| 1| 1| 1| 1
C% pst(*, 85): '2VJT6EA:    ' =>  0| 1| 1| 1| 1| 2|  |  
C% pst(*, 86): '2VJT6EB:    ' =>  0| 1| 1| 1| 1| 2|  |  
C% pst(*, 87): '2VJT6EC:    ' =>  0| 1| 1| 1| 1| 2|  |  
C% pst(*, 88): '2VEEAB:     ' =>  0| 1| 2| 4| 8|13|  |  
C% pst(*, 89): '2VEEAC:     ' =>  0| 1| 2| 4| 8|13|  |  
C% pst(*, 90): '2VEEBC:     ' =>  0| 1| 2| 4| 8|13|  |  
C% pst(*, 91): '2VEAA:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 92): '2VEBA:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 93): '2VECA:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 94): '2VEAB:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 95): '2VEBB:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 96): '2VECB:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 97): '2VEAC:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 98): '2VEBC:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*, 99): '2VECC:      ' =>  0| 0| 1| 2| 3| 4|  |  
C% pst(*,100): '2VAB:       ' =>  0| 1| 2| 3| 4| 5|  |  
C% pst(*,101): '2VAC:       ' =>  0| 1| 2| 3| 4| 5|  |  
C% pst(*,102): '2VBC:       ' =>  0| 1| 2| 3| 4| 5|  |  
C% pst(*,103): '2VAAA:      ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*,104): '2VEEE:      ' =>  0| 0| 1| 6|  |  |  |   
C% pst(*,105): '2VAEAEB:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*,106): '2VBEAEB:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*,107): '2VCEAEB:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*,108): '2VAEAEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*,109): '2VBEAEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*,110): '2VCEAEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*,111): '2VAEBEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*,112): '2VBEBEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*,113): '2VCEBEC:    ' =>  0| 0| 1| 3|  |  |  |  
C% pst(*,114): '2VABEA:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*,115): '2VACEA:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*,116): '2VBCEA:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*,117): '2VABEB:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*,118): '2VACEB:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*,119): '2VBCEB:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*,120): '2VABEC:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*,121): '2VACEC:     ' =>  0| 0| 0| 1|  |  |  |  
C% pst(*,122): '2VBCEC:     ' =>  0| 0| 0| 1|  |  |  |   
C% pst(*,123): '2VBLANCK1:  ' =>   |  |  |  |  |  |  |   
C% pst(*,124): '2VBLANCK2:  ' =>   |  |  |  |  |  |  |   
C% pst(*,125): '2VBLANCK3:  ' =>   |  |  |  |  |  |  |   
C% pst(*,126): '2VBLANCK4:  ' =>   |  |  |  |  |  |  |   
C% pst(*,127): '2VBLANCK5:  ' =>   |  |  |  |  |  |  |   
C% pst(*,128): '2VBLANCK6:  ' =>   |  |  |  |  |  |  |   
C% pst(*,129): '2VBLANCK7:  ' =>   |  |  |  |  |  |  |   
C% pst(*,130): '2VBLANCK8:  ' =>   |  |  |  |  |  |  |   
C% pst(*,131): '2VBLANCK9:  ' =>   |  |  |  |  |  |  |   
C% pst(*,132): '2VBLANCK10: ' =>   |  |  |  |  |  |  |   
C% pst(*,133): '2VBLANCK11: ' =>   |  |  |  |  |  |  |   
C% pst(*,134): '2VBLANCK12: ' =>   |  |  |  |  |  |  |   
C% pst(*,135): '2VBLANCK13: ' =>   |  |  |  |  |  |  |   
C% pst(*,136): '2VBLANCK14: ' =>   |  |  |  |  |  |  |   
C% pst(*,137): '2VBLANCK15: ' =>   |  |  |  |  |  |  |   
C% pst(*,138): '2VBLANCK16: ' =>   |  |  |  |  |  |  |   
C% pst(*,139): '2VBLANCK17: ' =>   |  |  |  |  |  |  |   
C% pst(*,140): '2VBLANCK18: ' =>   |  |  |  |  |  |  |   
C% pst(*,141): '2VBLANCK19: ' =>   |  |  |  |  |  |  |   
C% pst(*,142): '2VBLANCK20: ' =>   |  |  |  |  |  |  |   
C% pst(*,143): '2VBLANCK21: ' =>   |  |  |  |  |  |  |   
C% pst(*,144): '2VBLANCK22: ' =>   |  |  |  |  |  |  |   
C% pst(*,145): '2VBLANCK23: ' =>   |  |  |  |  |  |  |   
C% pst(*,146): '2VBLANCK24: ' =>   |  |  |  |  |  |  |   
C% pst(*,147): '2VBLANCK25: ' =>   |  |  |  |  |  |  |   
C% pst(*,148): '2VBLANCK26: ' =>   |  |  |  |  |  |  |   
C% pst(*,149): '2VBLANCK27: ' =>   |  |  |  |  |  |  |   
C% pst(*,150): '2VBLANCK28: ' =>   |  |  |  |  |  |  |   
C% pst(*,151): '2VBLANCK29: ' =>   |  |  |  |  |  |  |   
C% pst(*,152): '2VBLANCK30: ' =>   |  |  |  |  |  |  |  
C% pst(*,153): '1WJT6E4:    ' =>  1| 1| 1| 2| 2| 2|  |  
C% pst(*,154): '1WJT6E5:    ' =>  1| 1| 1| 2| 2| 2|  |  
C% pst(*,155): '1WJT6E6:    ' =>  1| 1| 1| 2| 2| 2|  |  
C% pst(*,156): '1WCOUPEEAB: ' =>  0| 1| 4| 9|12|23|  |  
C% pst(*,157): '1WCOUPEEAC: ' =>  0| 1| 4| 9|12|23|  |  
C% pst(*,158): '1WCOUPEEBC: ' =>  0| 1| 4| 9|12|23|  |  
C% pst(*,159): '1WCOUPEAA:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,160): '1WCOUPEBA:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,161): '1WCOUPECA:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,162): '1WCOUPEAB:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,163): '1WCOUPEBB:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,164): '1WCOUPECB:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,165): '1WCOUPEAC:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,166): '1WCOUPEBC:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,167): '1WCOUPECC:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,168): '1WCOUPAAA:  ' =>  0| 0| 0| 0|  |  |  |  
C% pst(*,169): '1WCOUPEEE:  ' =>  0| 0| 3|12|  |  |  |   
C% pst(*,170): '1WCOUAEAEB: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,171): '1WCOUBEAEB: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,172): '1WCOUCEAEB: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,173): '1WCOUAEAEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,174): '1WCOUBEAEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,175): '1WCOUCEAEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,176): '1WCOUAEBEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,177): '1WCOUBEBEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,178): '1WCOUCEBEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,179): '1WCOUPABEA: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,180): '1WCOUPACEA: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,181): '1WCOUPBCEA: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,182): '1WCOUPABEB: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,183): '1WCOUPACEB: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,184): '1WCOUPBCEB: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,185): '1WCOUPABEC: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,186): '1WCOUPACEC: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,187): '1WCOUPBCEC: ' =>  0| 0| 1| 1|  |  |  |   
C% pst(*,188): '1WBLANCK1:  ' =>   |  |  |  |  |  |  |   
C% pst(*,189): '1WBLANCK2:  ' =>   |  |  |  |  |  |  |   
C% pst(*,190): '1WBLANCK3:  ' =>   |  |  |  |  |  |  |   
C% pst(*,191): '1WBLANCK4:  ' =>   |  |  |  |  |  |  |   
C% pst(*,192): '1WBLANCK5:  ' =>   |  |  |  |  |  |  |   
C% pst(*,193): '1WBLANCK6:  ' =>   |  |  |  |  |  |  |   
C% pst(*,194): '1WBLANCK7:  ' =>   |  |  |  |  |  |  |   
C% pst(*,195): '1WBLANCK8:  ' =>   |  |  |  |  |  |  |   
C% pst(*,196): '1WBLANCK9:  ' =>   |  |  |  |  |  |  |   
C% pst(*,197): '1WBLANCK10: ' =>   |  |  |  |  |  |  |   
C% pst(*,198): '1WBLANCK11: ' =>   |  |  |  |  |  |  |   
C% pst(*,199): '1WBLANCK12: ' =>   |  |  |  |  |  |  |   
C% pst(*,200): '1WBLANCK13: ' =>   |  |  |  |  |  |  |   
C% pst(*,201): '1WBLANCK14: ' =>   |  |  |  |  |  |  |   
C% pst(*,202): '1WBLANCK15: ' =>   |  |  |  |  |  |  |   
C% pst(*,203): '1WBLANCK16: ' =>   |  |  |  |  |  |  |   
C% pst(*,204): '1WBLANCK17: ' =>   |  |  |  |  |  |  |   
C% pst(*,205): '1WBLANCK18: ' =>   |  |  |  |  |  |  |   
C% pst(*,206): '1WBLANCK19: ' =>   |  |  |  |  |  |  |   
C% pst(*,207): '1WBLANCK20: ' =>   |  |  |  |  |  |  |   
C% pst(*,208): '1WBLANCK21: ' =>   |  |  |  |  |  |  |   
C% pst(*,209): '1WBLANCK22: ' =>   |  |  |  |  |  |  |   
C% pst(*,210): '1WBLANCK23: ' =>   |  |  |  |  |  |  |   
C% pst(*,211): '1WBLANCK24: ' =>   |  |  |  |  |  |  |   
C% pst(*,212): '1WBLANCK25: ' =>   |  |  |  |  |  |  |   
C% pst(*,213): '1WBLANCK26: ' =>   |  |  |  |  |  |  |   
C% pst(*,214): '1WBLANCK27: ' =>   |  |  |  |  |  |  |   
C% pst(*,215): '1WBLANCK28: ' =>   |  |  |  |  |  |  |   
C% pst(*,216): '1WBLANCK29: ' =>   |  |  |  |  |  |  |   
C% pst(*,217): '1WBLANCK30: ' =>   |  |  |  |  |  |  |  
C% pst(*,218): '2WJT6E4:    ' =>  1| 1| 1| 2| 2| 2|  |  
C% pst(*,219): '2WJT6E5:    ' =>  1| 1| 1| 2| 2| 2|  |  
C% pst(*,220): '2WJT6E6:    ' =>  1| 1| 1| 2| 2| 2|  |  
C% pst(*,221): '2WCOUPEEAB: ' =>  0| 1| 4| 9|12|23|  |  
C% pst(*,222): '2WCOUPEEAC: ' =>  0| 1| 4| 9|12|23|  |  
C% pst(*,223): '2WCOUPEEBC: ' =>  0| 1| 4| 9|12|23|  |  
C% pst(*,224): '2WCOUPEAA:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,225): '2WCOUPEBA:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,226): '2WCOUPECA:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,227): '2WCOUPEAB:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,228): '2WCOUPEBB:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,229): '2WCOUPECB:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,230): '2WCOUPEAC:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,231): '2WCOUPEBC:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,232): '2WCOUPECC:  ' =>  0| 1| 2| 3| 5| 7|  |  
C% pst(*,233): '2WCOUPAAA:  ' =>  0| 0| 0| 0|  |  |  |  
C% pst(*,234): '2WCOUPEEE:  ' =>  0| 0| 3|12|  |  |  |   
C% pst(*,235): '2WCOUAEAEB: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,236): '2WCOUBEAEB: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,237): '2WCOUCEAEB: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,238): '2WCOUAEAEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,239): '2WCOUBEAEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,240): '2WCOUCEAEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,241): '2WCOUAEBEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,242): '2WCOUBEBEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,243): '2WCOUCEBEC: ' =>  0| 0| 1| 5|  |  |  |  
C% pst(*,244): '2WCOUPABEA: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,245): '2WCOUPACEA: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,246): '2WCOUPBCEA: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,247): '2WCOUPABEB: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,248): '2WCOUPACEB: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,249): '2WCOUPBCEB: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,250): '2WCOUPABEC: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,251): '2WCOUPACEC: ' =>  0| 0| 1| 1|  |  |  |  
C% pst(*,252): '2WCOUPBCEC: ' =>  0| 0| 1| 1|  |  |  |    
C% pst(*,253): '2WBLANCK1:  ' =>   |  |  |  |  |  |  |   
C% pst(*,254): '2WBLANCK2:  ' =>   |  |  |  |  |  |  |   
C% pst(*,255): '2WBLANCK3:  ' =>   |  |  |  |  |  |  |   
C% pst(*,256): '2WBLANCK4:  ' =>   |  |  |  |  |  |  |   
C% pst(*,257): '2WBLANCK5:  ' =>   |  |  |  |  |  |  |   
C% pst(*,258): '2WBLANCK6:  ' =>   |  |  |  |  |  |  |   
C% pst(*,259): '2WBLANCK7:  ' =>   |  |  |  |  |  |  |   
C% pst(*,260): '2WBLANCK8:  ' =>   |  |  |  |  |  |  |   
C% pst(*,261): '2WBLANCK9:  ' =>   |  |  |  |  |  |  |   
C% pst(*,262): '2WBLANCK10: ' =>   |  |  |  |  |  |  |   
C% pst(*,263): '2WBLANCK11: ' =>   |  |  |  |  |  |  |   
C% pst(*,264): '2WBLANCK12: ' =>   |  |  |  |  |  |  |   
C% pst(*,265): '2WBLANCK13: ' =>   |  |  |  |  |  |  |   
C% pst(*,266): '2WBLANCK14: ' =>   |  |  |  |  |  |  |   
C% pst(*,267): '2WBLANCK15: ' =>   |  |  |  |  |  |  |   
C% pst(*,268): '2WBLANCK16: ' =>   |  |  |  |  |  |  |   
C% pst(*,269): '2WBLANCK17: ' =>   |  |  |  |  |  |  |   
C% pst(*,270): '2WBLANCK18: ' =>   |  |  |  |  |  |  |   
C% pst(*,271): '2WBLANCK19: ' =>   |  |  |  |  |  |  |   
C% pst(*,272): '2WBLANCK20: ' =>   |  |  |  |  |  |  |   
C% pst(*,273): '2WBLANCK21: ' =>   |  |  |  |  |  |  |   
C% pst(*,274): '2WBLANCK22: ' =>   |  |  |  |  |  |  |   
C% pst(*,275): '2WBLANCK23: ' =>   |  |  |  |  |  |  |   
C% pst(*,276): '2WBLANCK24: ' =>   |  |  |  |  |  |  |   
C% pst(*,277): '2WBLANCK25: ' =>   |  |  |  |  |  |  |   
C% pst(*,278): '2WBLANCK26: ' =>   |  |  |  |  |  |  |   
C% pst(*,279): '2WBLANCK27: ' =>   |  |  |  |  |  |  |   
C% pst(*,280): '2WBLANCK28: ' =>   |  |  |  |  |  |  |   
C% pst(*,281): '2WBLANCK29: ' =>   |  |  |  |  |  |  |   
C% pst(*,282): '2WBLANCK30: ' =>   |  |  |  |  |  |  |   
C% pst(*,283): 'SO_O0:      ' =>  2|  |  |  |  |  |  |   
C% pst(*,284): 'SO_O1_A:    ' =>   | 5|  |  |  |  |  |   
C% pst(*,285): 'SO_O2_A:    ' =>   |  |10|  |  |  |  |   
C% pst(*,286): 'SO_O3_A:    ' =>   |  |  |16|  |  |  |   
C% pst(*,287): 'SO_O4_A:    ' =>   |  |  |  |24|  |  |   
C% pst(*,288): 'SO_O1_B:    ' =>   | 5|  |  |  |  |  |  
C% pst(*,289): 'SO_O2_B:    ' =>   |  |10|  |  |  |  |  
C% pst(*,290): 'SO_O3_B:    ' =>   |  |  |16|  |  |  |  
C% pst(*,291): 'SO_O4_B:    ' =>   |  |  |  |24|  |  |  
C% pst(*,292): 'SO_O1_C:    ' =>   | 5|  |  |  |  |  |  
C% pst(*,293): 'SO_O2_C:    ' =>   |  |10|  |  |  |  |  
C% pst(*,294): 'SO_O3_C:    ' =>   |  |  |16|  |  |  |  
C% pst(*,295): 'SO_O4_C:    ' =>   |  |  |  |24|  |  |  
C% pst(*,296): 'SO_O2_AB:   ' =>   |  |15|  |  |  |  |   
C% pst(*,297): 'SO_O3_AB:   ' =>   |  |  |48|  |  |  |   
C% pst(*,298): 'SO_O4_AB:   ' =>   |  |  |  |  |  |  |   
C% pst(*,299): 'SO_O2_AC:   ' =>   |  |15|  |  |  |  |   
C% pst(*,300): 'SO_O3_AC:   ' =>   |  |  |48|  |  |  |   
C% pst(*,301): 'SO_O4_AC:   ' =>   |  |  |  |  |  |  |   
C% pst(*,302): 'SO_O2_BC:   ' =>   |  |15|  |  |  |  |   
C% pst(*,303): 'SO_O3_BC:   ' =>   |  |  |48|  |  |  |   
C% pst(*,304): 'SO_O4_BC:   ' =>   |  |  |  |  |  |  |  
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------
! nstat:	dimension of matrix
! q:	coordinate
! p:	parameter
! e:	diabatic energy matrix

      subroutine diab(e,n,q,p,pst,iref,npar)
      use omp_lib
      implicit none

      include 'params.incl'
      include 'states.incl'

      integer npar !number of parameters
      integer n, i, j, pst(2,np), iref
      double precision e(nstat,nstat), q(qn), p(npar)
      double precision w
      double precision lq(qn)
      logical dbg
      double precision r(qn) !for ctrans
      logical skip           !for ctrans

      dbg = .true.
      dbg = .false.

!..   compute potential elements every time:
      j=omp_get_thread_num()

      do i=1,9
        lq(i)=q(i)
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! killed without tmcs !!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Transformation of coordinates
      call ctrans(q,lq,r,1,skip,p,npar)

c     Berechnung der Monome
      call vwzprec(n,lq)

      !nullen der diabatischen Matrix
      do i=1, nstat
         do j=1, nstat
           e(i,j) = 0.d0
         enddo
      enddo

      !---------------------------------------------------------------
      !UNCOUPLED
      !Calculation of the uncoupled diagonalelements of the Matrix
      !without JT and PJT
      !
      !Only the relevent part of the pst-Vektor is given, to make the
      !function more generic
      !
      !Exitation energy is given
      !
      !Builds up_ the a-state element
      call vstat(e(3,3),n,p,pst(1,11),iref,p(pst(1,1)+1),pst(1,81))
      !builds up_ the e-state element
      call vstat(e(1,1),n,p,pst(1,82),iref,p(pst(1,1)),pst(1,152))
      e(2,2) = e(1,1)

      !---------------------------------------------------------------
      !JT-COUPLING
      !
      !Only the relevent part of the pst-Vektor is given, to make the
      !function more generic
      !
      !Builds up_ the JT-coupling for the diagonal e x e elements
      call wcoup(w,n,p,pst(1,153),iref,pst(1,217))
      e(1,1) = e(1,1) - w
      e(2,2) = e(2,2) + w
      !Builds up_ the JT-coupling for the ofdiagonal e x e elements
      call zcoup(e(1,2),n,p,pst(1,153),iref,pst(1,217))
      e(2,1) = e(1,2)
      
      !---------------------------------------------------------------
      !PJT-COUPLING
      !
      !Only the relevent part of the pst-Vektor is given, to make the
      !function more generic
      !
      !Builds up_ the PJT-coupling
      call wcoup(e(3,1),n,p,pst(1,218),iref,pst(1,282))
      e(1,3) = e(3,1)
      call zcoup(e(2,3),n,p,pst(1,218),iref,pst(1,282))
      e(3,2) = e(2,3)

100   format(5f14.7)
      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE VSTAT(E,N,P,PST,IREF,VEXEN,NPAR)
C%
C% subroutine to create E state matrix element this corresponds to 
C% morse for sym stretch, and V of JT6 for e modes
C%
C% you should give only the part of the pst which starts with the 
C% used block markers
C% pst(*, 1): epolynoma 
C% pst(*, 2): epolynomb 
C% pst(*, 3): epolynomc 
C% pst(*, 4): evjt6a    
C% pst(*, 5): evjt6b    
C% pst(*, 6): evjt6c    
C% pst(*, 7): ev_eeab   
C% pst(*, 8): ev_eeac 
C% pst(*, 9): ev_eebc   
C% pst(*,10): ev_eaa    
C% pst(*,11): ev_eba    
C% pst(*,12): ev_eca    
C% pst(*,13): ev_eab    
C% pst(*,14): ev_ebb    
C% pst(*,15): ev_ecb    
C% pst(*,16): ev_eac    
C% pst(*,17): ev_ebc    
C% pst(*,18): ev_ecc    
C% pst(*,19): ev_aa     
C% pst(*,20): ev_aa     
C% pst(*,21): ev_aa
C%
C% Input variables:
C% n:     number of the current point (int)
C% p:     parametervector (double[npar])
C% pst:   pointers for parameters
C% iref:  1: ref. model, 0 correction (int)
C% vexen: exitation energy (double)
C% npar:  number of given parameters (int)
C% 
C% Output variables:
C% e: matrix elements (double)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine vstat(e,n,p,pst,iref,vexen,npar)
      implicit none

      include 'params.incl'
      include 'states.incl'

      integer n !runing indices

      integer pst(2,nv)     !position_ of parameters
      integer npar          !number of given parameters
      double precision p(npar) !Parameters 

      integer iref !if it is the refference Hamiltonien

      double precision e !energy

      double precision dum !Matrix element

      double precision vexen !Vertical excitation energy

      !write(6,'(''In vstat'')')
      !write(6,'(5f9.5)') p 

!.....vertical excitation energy:
      dum=0.d0
      e = vexen
      if (iref.ne.0) e = 0.d0

!..   single mode uncoupled diagonal potentials:
      !Buidlds up_ the uncoupled a-modes:
      call epolynoma(dum,n, p(pst(1,1)), pst(2,1), iref)
      e=e+dum
      call epolynomb(dum,n, p(pst(1,2)), pst(2,2), iref)
      e=e+dum
      call epolynomc(dum,n, p(pst(1,3)), pst(2,3), iref)
      e=e+dum

      !Builds up_ the uncoupled e-modes
      call evjt6a(dum,n, p(pst(1,4)), pst(2,4), iref)
      e = e + dum
      call evjt6b(dum,n, p(pst(1,5)), pst(2,5), iref) 
      e = e + dum
      call evjt6c(dum,n, p(pst(1,6)), pst(2,6), iref) 
      e = e + dum

!..   mode-mode coupling e-e
      call ev_eeab(dum,n, p(pst(1,7)), pst(2,7), iref)
      e = e + dum
      call ev_eeac(dum,n, p(pst(1,8)), pst(2,8), iref)
      e = e + dum
      call ev_eebc(dum,n, p(pst(1,9)), pst(2,9), iref)
      e = e + dum

!..   mode-mode coupling of a1-e
      call ev_eaa(dum,n, p(pst(1,10)), pst(2,10), iref)  
      e = e + dum                       
      call ev_eba(dum,n, p(pst(1,11)), pst(2,11), iref)  
      e = e + dum                              
      call ev_eca(dum,n, p(pst(1,12)), pst(2,12), iref)  
      e = e + dum                              
      call ev_eab(dum,n, p(pst(1,13)), pst(2,13), iref)  
      e = e + dum                              
      call ev_ebb(dum,n, p(pst(1,14)), pst(2,14), iref)  
      e = e + dum                              
      call ev_ecb(dum,n, p(pst(1,15)), pst(2,15), iref)  
      e = e + dum
      call ev_eac(dum,n, p(pst(1,16)), pst(2,16), iref)  
      e = e + dum                              
      call ev_ebc(dum,n, p(pst(1,17)), pst(2,17), iref)  
      e = e + dum                              
      call ev_ecc(dum,n, p(pst(1,18)), pst(2,18), iref)  
      e = e + dum

!..   mode-mode coupling of a1-a1
      call Ev_ab(dum,n,p(pst(1,19)),pst(2,19), iref)
      e = e + dum
      call Ev_ac(dum,n,p(pst(1,20)),pst(2,20), iref)
      e = e + dum
      call Ev_bc(dum,n,p(pst(1,21)),pst(2,21), iref)
      e = e + dum

c     3-mode-coupling
      call e_vaaa(dum,n,p(pst(1,22)),pst(2,22),iref)
      e = e + dum
      call e_veee(dum,n,p(pst(1,23)),pst(2,23),iref)
      e = e + dum

      call e_vaeaeb(dum,n,p(pst(1,24)),pst(2,24),iref)
      e = e + dum
      call e_vbeaeb(dum,n,p(pst(1,25)),pst(2,25),iref)
      e = e + dum
      call e_vceaeb(dum,n,p(pst(1,26)),pst(2,26),iref)
      e = e + dum
      call e_vaeaec(dum,n,p(pst(1,27)),pst(2,27),iref)
      e = e + dum
      call e_vbeaec(dum,n,p(pst(1,28)),pst(2,28),iref)
      e = e + dum
      call e_vceaec(dum,n,p(pst(1,29)),pst(2,29),iref)
      e = e + dum
      call e_vaebec(dum,n,p(pst(1,30)),pst(2,30),iref)
      e = e + dum
      call e_vbebec(dum,n,p(pst(1,31)),pst(2,31),iref)
      e = e + dum
      call e_vcebec(dum,n,p(pst(1,32)),pst(2,32),iref)
      e = e + dum

      call e_vabea(dum,n,p(pst(1,33)),pst(2,33),iref)
      e = e + dum
      call e_vacea(dum,n,p(pst(1,34)),pst(2,34),iref)
      e = e + dum
      call e_vbcea(dum,n,p(pst(1,35)),pst(2,35),iref)
      e = e + dum
      call e_vabeb(dum,n,p(pst(1,36)),pst(2,36),iref)
      e = e + dum
      call e_vaceb(dum,n,p(pst(1,37)),pst(2,37),iref)
      e = e + dum
      call e_vbceb(dum,n,p(pst(1,38)),pst(2,38),iref)
      e = e + dum
      call e_vabec(dum,n,p(pst(1,39)),pst(2,39),iref)
      e = e + dum
      call e_vacec(dum,n,p(pst(1,40)),pst(2,40),iref)
      e = e + dum
      call e_vbcec(dum,n,p(pst(1,41)),pst(2,41),iref)
      e = e + dum

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c subroutine to create W matrix elements for JT and PJT
c
c you should give only the part of the pst which starts with theused
c block markers
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoup(e,n,p,pst,iref,npar)
      implicit none

      include 'params.incl'
      include 'states.incl'

      integer n !runing indices

      integer pst(2,nwz)       !position_ of parameters
      integer npar             !number of parameters
      double precision p(npar) !Parameters 

      integer iref !if it is the refference Hamiltonien

      double precision e !energy

      double precision dum !Matrix element

!.....vertical excitation energy:
      e = 0.d0
      if (iref.eq.1) e = 0.d0

!..   add (P)JT coupling for 1. e mode
      call Ewjt6a(dum,n, p(pst(1,1)), pst(2,1), iref)
      e=e+dum
      call Ewjt6b(dum,n, p(pst(1,2)), pst(2,2), iref)
      e=e+dum
      call Ewjt6c(dum,n, p(pst(1,3)), pst(2,3), iref)
      e=e+dum
 
!..   add (P)JT mode-mode coupling between two e modes:
      call Ewcoup_eeab(dum,n,p(pst(1,4)), pst(2,4), iref)
      e=e+dum
      call Ewcoup_eeac(dum,n,p(pst(1,5)), pst(2,5), iref)
      e=e+dum
      call Ewcoup_eebc(dum,n,p(pst(1,6)), pst(2,6), iref)
      e=e+dum
 
!..   add (P)JT mode-mode coupling between a1 and e modes:
      call Ewcoup_eaa(dum,n, p(pst(1, 7)), pst(2, 7), iref)
      e=e+dum                                      
      call Ewcoup_eba(dum,n, p(pst(1, 8)), pst(2, 8), iref)
      e=e+dum                                      
      call Ewcoup_eca(dum,n, p(pst(1, 9)), pst(2, 9), iref)
      e=e+dum                                      
      call Ewcoup_eab(dum,n, p(pst(1,10)), pst(2,10), iref)
      e=e+dum                                      
      call Ewcoup_ebb(dum,n, p(pst(1,11)), pst(2,11), iref)
      e=e+dum                                      
      call Ewcoup_ecb(dum,n, p(pst(1,12)), pst(2,12), iref)
      e=e+dum                                      
      call Ewcoup_eac(dum,n, p(pst(1,13)), pst(2,13), iref)
      e=e+dum                                      
      call Ewcoup_ebc(dum,n, p(pst(1,14)), pst(2,14), iref)
      e=e+dum                                      
      call Ewcoup_ecc(dum,n, p(pst(1,15)), pst(2,15), iref)
      e=e+dum

c     3-mode-coupling
      call wcoupeee(dum,n,p(pst(1,17)),pst(2,17),iref)
      e=e+dum

      call wcoupaeaeb(dum,n,p(pst(1,18)),pst(2,18),iref)
      e=e+dum
      call wcoupbeaeb(dum,n,p(pst(1,19)),pst(2,19),iref)
      e=e+dum
      call wcoupceaeb(dum,n,p(pst(1,20)),pst(2,20),iref)
      e=e+dum
      call wcoupaeaec(dum,n,p(pst(1,21)),pst(2,21),iref)
      e=e+dum
      call wcoupbeaec(dum,n,p(pst(1,22)),pst(2,22),iref)
      e=e+dum
      call wcoupceaec(dum,n,p(pst(1,23)),pst(2,23),iref)
      e=e+dum
      call wcoupaebec(dum,n,p(pst(1,24)),pst(2,24),iref)
      e=e+dum
      call wcoupbebec(dum,n,p(pst(1,25)),pst(2,25),iref)
      e=e+dum
      call wcoupcebec(dum,n,p(pst(1,26)),pst(2,26),iref)
      e=e+dum

      call wcoupabea(dum,n,p(pst(1,27)),pst(2,27),iref)
      e=e+dum
      call wcoupacea(dum,n,p(pst(1,28)),pst(2,28),iref)
      e=e+dum
      call wcoupbcea(dum,n,p(pst(1,29)),pst(2,29),iref)
      e=e+dum
      call wcoupabeb(dum,n,p(pst(1,30)),pst(2,30),iref)
      e=e+dum
      call wcoupaceb(dum,n,p(pst(1,31)),pst(2,31),iref)
      e=e+dum
      call wcoupbceb(dum,n,p(pst(1,32)),pst(2,32),iref)
      e=e+dum
      call wcoupabec(dum,n,p(pst(1,33)),pst(2,33),iref)
      e=e+dum
      call wcoupacec(dum,n,p(pst(1,34)),pst(2,34),iref)
      e=e+dum
      call wcoupbcec(dum,n,p(pst(1,35)),pst(2,35),iref)
      e=e+dum

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c subroutine to create Z matrix elements for JT and PJT
c
c you should give only the part of the pst which starts with theused
c block markers
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoup(e,n,p,pst,iref,npar)
      implicit none

      include 'params.incl'
      include 'states.incl'

      integer n !runing indices

      integer pst(2,nwz)       !position_ of parameters
      integer npar             !number of parameters
      double precision p(npar) !Parameters 

      integer iref !if it is the refference Hamiltonien

      double precision e !energy

      double precision dum !Matrix element

c     vertical excitation energy:
      e = 0.d0
      if (iref.eq.1) e = 0.d0

c     add (P)JT coupling for 1. e mode
      call Ezjt6a(dum,n, p(pst(1,1)), pst(2,1), iref)
      e = e+dum
      call Ezjt6b(dum,n, p(pst(1,2)), pst(2,2), iref)
      e = e+dum
      call Ezjt6c(dum,n, p(pst(1,3)), pst(2,3), iref)
      e = e+dum

c     add (P)JT mode-mode coupling between two e modes:
      call Ezcoup_eeab(dum,n,p(pst(1,4)), pst(2,4), iref)
      e=e+dum
      call Ezcoup_eeac(dum,n,p(pst(1,5)), pst(2,5), iref)
      e=e+dum
      call Ezcoup_eebc(dum,n,p(pst(1,6)), pst(2,6), iref)
      e=e+dum

c     add (P)JT mode-mode coupling between a1 and e modes:
      call Ezcoup_eaa(dum,n, p(pst(1, 7)), pst(2, 7), iref)
      e=e+dum                                      
      call Ezcoup_eba(dum,n, p(pst(1, 8)), pst(2, 8), iref)
      e=e+dum                                      
      call Ezcoup_eca(dum,n, p(pst(1, 9)), pst(2, 9), iref)
      e=e+dum                                      
      call Ezcoup_eab(dum,n, p(pst(1,10)), pst(2,10), iref)
      e=e+dum                                      
      call Ezcoup_ebb(dum,n, p(pst(1,11)), pst(2,11), iref)
      e=e+dum                                      
      call Ezcoup_ecb(dum,n, p(pst(1,12)), pst(2,12), iref)
      e=e+dum                                      
      call Ezcoup_eac(dum,n, p(pst(1,13)), pst(2,13), iref)
      e=e+dum                                      
      call Ezcoup_ebc(dum,n, p(pst(1,14)), pst(2,14), iref)
      e=e+dum                                      
      call Ezcoup_ecc(dum,n, p(pst(1,15)), pst(2,15), iref)
      e=e+dum

c     3-mode-coupling
      call zcoupeee(dum,n,p(pst(1,17)),pst(2,17),iref)
      e=e+dum

      call zcoupaeaeb(dum,n,p(pst(1,18)),pst(2,18),iref)
      e=e+dum
      call zcoupbeaeb(dum,n,p(pst(1,19)),pst(2,19),iref)
      e=e+dum
      call zcoupceaeb(dum,n,p(pst(1,20)),pst(2,20),iref)
      e=e+dum
      call zcoupaeaec(dum,n,p(pst(1,21)),pst(2,21),iref)
      e=e+dum
      call zcoupbeaec(dum,n,p(pst(1,22)),pst(2,22),iref)
      e=e+dum
      call zcoupceaec(dum,n,p(pst(1,23)),pst(2,23),iref)
      e=e+dum
      call zcoupaebec(dum,n,p(pst(1,24)),pst(2,24),iref)
      e=e+dum
      call zcoupbebec(dum,n,p(pst(1,25)),pst(2,25),iref)
      e=e+dum
      call zcoupcebec(dum,n,p(pst(1,26)),pst(2,26),iref)
      e=e+dum

      call zcoupabea(dum,n,p(pst(1,27)),pst(2,27),iref)
      e=e+dum
      call zcoupacea(dum,n,p(pst(1,28)),pst(2,28),iref)
      e=e+dum
      call zcoupbcea(dum,n,p(pst(1,29)),pst(2,29),iref)
      e=e+dum
      call zcoupabeb(dum,n,p(pst(1,30)),pst(2,30),iref)
      e=e+dum
      call zcoupaceb(dum,n,p(pst(1,31)),pst(2,31),iref)
      e=e+dum
      call zcoupbceb(dum,n,p(pst(1,32)),pst(2,32),iref)
      e=e+dum
      call zcoupabec(dum,n,p(pst(1,33)),pst(2,33),iref)
      e=e+dum
      call zcoupacec(dum,n,p(pst(1,34)),pst(2,34),iref)
      e=e+dum
      call zcoupbcec(dum,n,p(pst(1,35)),pst(2,35),iref)
      e=e+dum

      end
       
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Builds up_ the a-Mode up_ to 8th order
c uses not more than 8 parameters
c 1.: 1
c 2.: 1
c 3.: 1
c 4.: 1
c 5.: 1
c 6.: 1
c 7.: 1
c 8.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine epolynoma(polynoma,n, p, j, iref)
      implicit none

      integer i !running indicies

      integer n !point
      
      integer j               !number of given parameters
      double precision p(j)   !internal parameter vector


      integer ii, jj, iref !used for ref. Hamiltonien

      double precision polynoma !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

      if (j.gt.8) then
        write(6,*) 'Error: polynom higher than 8th order not implemente
     &d!'
        stop
      endif
      
      polynoma = 0.d0

      ii=1
      jj=min(2,j)
      if (iref.eq.1) then
        ii=3
        jj=j
      endif

      do i = ii, jj
         polynoma = polynoma + pa(i,n) * p(i)
      enddo

      end
       
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Builds up_ the b-Mode up_ to 8th order
c uses not more than 8 parameters
c 1.: 1
c 2.: 1
c 3.: 1
c 4.: 1
c 5.: 1
c 6.: 1
c 7.: 1
c 8.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine epolynomb(polynoma,n, p, j, iref)
      implicit none

      integer i !running indicies

      integer n !point
      
      integer j               !number of given parameters
      double precision p(j) !internal parameter vector


      integer ii, jj, iref !used for ref. Hamiltonien

      double precision polynoma !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

      if (j.gt.8) then
        write(6,*) 'Error: polynom higher than 8th order not implemente
     &d!'
        stop
      endif
      
      polynoma = 0.d0

      ii=1
      jj=min(2,j)
      if (iref.eq.1) then
        ii=3
        jj=j
      endif

      do i = ii, jj
         polynoma = polynoma + pb(i,n) * p(i) 
      enddo

      end
       
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Builds up_ the c-Mode up_ to 8th order
c uses not more than 8 parameters
c 1.: 1
c 2.: 1
c 3.: 1
c 4.: 1
c 5.: 1
c 6.: 1
c 7.: 1
c 8.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine epolynomc(polynoma,n, p, j, iref)
      implicit none

      integer i !running indicies

      integer n !point
      
      integer j               !number of given parameters
      double precision p(j) !internal parameter vector


      integer ii, jj, iref !used for ref. Hamiltonien

      double precision polynoma !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

      if (j.gt.8) then
        write(6,*) 'Error: polynom higher than 8th order not implemente
     &d!'
        stop
      endif
      
      polynoma = 0.d0

      ii=1
      jj=min(2,j)
      if (iref.eq.1) then
        ii=3
        jj=j
      endif

      do i = ii, jj
         polynoma = polynoma + pc(i,n) * p(i) 
      enddo

      end


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Builds up_ the first e-Mode up_ to sixth order
c uses not more than 6 parameters
c 2.: 1
c 3.: 1
c 4.: 1
c 5.: 1
C 6.: 2
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Evjt6a(vjt6a,n, par_, j, iref)
      implicit none

      integer i !running indicies

      integer n !point

      integer j               !number of given parameters
      integer lnx             !maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii, jj, iref !used for ref. Hamiltonien

      double precision vjt6a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1 
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo
      

!     2nd order
      vjt6a = p(1) * va(1,n) 

!     3rd order
      vjt6a = vjt6a + (p(2) * va(2,n))    

!     4th order
      vjt6a = vjt6a + (p(3) * va(3,n)) 

!     5th order
      vjt6a = vjt6a + (p(4) * va(4,n)) 

!     6th order
      vjt6a = vjt6a + (p(5) * va(5,n) + p(6) * va(6,n)) 
      end
     
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Builds up_ the second e-Mode up_ to sixth order
c uses not more than 6 parameters
c 2.: 1
c 3.: 1
c 4.: 1
c 5.: 1
C 6.: 2
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Evjt6b(vjt6b, n, par_, j, iref)
      implicit none

      integer i !running indicies

      integer n !point
      
      integer j               !number of given parameters
      integer lnx             !maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii, jj, iref !used for ref. Hamiltonien

      double precision vjt6b !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo
c     this takes care of split reference plus correction scheme
      ii=1 
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo

c     2nd order
      vjt6b = p(1) * vb(1,n)   

c     3rd order
      vjt6b = vjt6b + (p(2) * vb(2,n))     

c     4th order
      vjt6b = vjt6b + (p(3) * vb(3,n))    

c     5th order
      vjt6b = vjt6b + (p(4) * vb(4,n))  

c     6th order
      vjt6b = vjt6b + (p(5) * vb(5,n) + p(6) * vb(6,n))
     
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Builds up_ the third e-Mode up_ to sixth order
c uses not more than 6 parameters
c 2.: 1
c 3.: 1
c 4.: 1
c 5.: 1
C 6.: 2
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Evjt6c(vjt6c, n, par_, j, iref)
      implicit none

      integer i !running indicies
      
      integer n !point
      
      integer j               !number of given parameters
      integer lnx             !maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii, jj, iref !used for ref. Hamiltonien

      double precision vjt6c !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo
c     this takes care of split reference plus correction scheme
      ii=1 
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo

!     2nd order
      vjt6c = p(1) * vc(1,n)   

!     3rd order
      vjt6c = vjt6c + (p(2) * vc(2,n))    

!     4th order
      vjt6c = vjt6c + (p(3) * vc(3,n))   

!     5th order
      vjt6c = vjt6c + (p(4) * vc(4,n)) 

!     6th order
      vjt6c = vjt6c + (p(5) * vc(5,n) + p(6) * vc(6,n)) 
     
      end
     
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W Jahn-Teller matrix elements up_ to 6th order
c for ea
c uses not more then 9 Parameters
c 1.: 1
c 2.: 1
c 3.: 1
c 4.: 2
c 5.: 2
c 6.: 2
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewjt6a(wjt6a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running index

      integer j !Number of given Parameters
      integer lnx !Maximum number of Parameters
      parameter (lnx=9)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector
 

      integer iref, ii, jj !used for ref. Hamiltonien

      double precision wjt6a !Matrix elemenet

      include 'params.incl'
      include 'vwz_MeX.incl'

!      write(6,*), wa(1:9,n)
c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo

!     1st order 
      wjt6a = p(1) * wa(1,n)

!     2nd order 
      wjt6a = wjt6a + (p(2) * wa(2,n)) 

!     3rd order 
      wjt6a = wjt6a + (p(3) * wa(3,n)) 

!     4th order 
      wjt6a = wjt6a + (p(4) * wa(4,n) + p(5) * wa(5,n))

!     5th order 
      wjt6a = wjt6a + (p(6) * wa(6,n) + p(7) * wa(7,n))

!     6th order 
      wjt6a = wjt6a + (p(8) * wa(8,n) + p(9)* wa(9,n)) 

      end     

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W Jahn-Teller matrix elements up_ to 6th order
c for eb 
c uses not more then 9 Parameters
c 1.: 1
c 2.: 1
c 3.: 1
c 4.: 2
c 5.: 2
c 6.: 2
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewjt6b(wjt6a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running index

      integer j !Number of given Parameters
      integer lnx !Maximum number of Parameters
      parameter (lnx=9)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector
 

      integer iref, ii, jj !used for ref. Hamiltonien

      double precision wjt6a !Matrix elemenet

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo

!     1st order 
      wjt6a = p(1) * wb(1,n)

!     2nd order 
      wjt6a = wjt6a + (p(2) * wb(2,n)) 

!     3rd order 
      wjt6a = wjt6a + (p(3) * wb(3,n)) 

!     4th order 
      wjt6a = wjt6a + (p(4) * wb(4,n) + p(5) * wb(5,n)) 

!     5th order 
      wjt6a = wjt6a + (p(6) * wb(6,n) + p(7) * wb(7,n)) 

!     6th order 
      wjt6a = wjt6a + (p(8) * wb(8,n) + p(9)* wb(9,n))  

      end
     
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W Jahn-Teller matrix elements up_ to 6th order
c for ec
c uses not more then 9 Parameters
c 1.: 1
c 2.: 1
c 3.: 1
c 4.: 2
c 5.: 2
c 6.: 2
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewjt6c(wjt6a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running index

      integer j !Number of given Parameters
      integer lnx !Maximum number of Parameters
      parameter (lnx=9)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector
 

      integer iref, ii, jj !used for ref. Hamiltonien

      double precision wjt6a !Matrix elemenet

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo

!     1st order 
      wjt6a = p(1) * wc(1,n)

!     2nd order 
      wjt6a = wjt6a + (p(2) * wc(2,n)) 

!     3rd order 
      wjt6a = wjt6a + (p(3) * wc(3,n)) 

!     4th order 
      wjt6a = wjt6a + (p(4) * wc(4,n) + p(5) * wc(5,n)) 

!     5th order 
      wjt6a = wjt6a + (p(6) * wc(6,n) + p(7) * wc(7,n)) 

!     6th order 
      wjt6a = wjt6a + (p(8) * wc(8,n) + p(9)* wc(9,n)) 

      end
 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z Jahn-Teller matrix elements up_ to 6th order
c for ea
c uses not more then 9 Parameters
c 1.: 1
c 2.: 1
c 3.: 1
c 4.: 2
c 5.: 2
c 6.: 2
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezjt6a(zjt6a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running index

      integer j               !Number of given Parameters
      integer lnx             !Maximum number of Parameters
      parameter (lnx=9)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector
 

      integer iref, ii, jj !used for ref. Hamiltonien

      double precision zjt6a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo
      
c     this takes care of split reference plus correction scheme
      ii=1
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo

!     1st order
      zjt6a = p(1) * za(1,n)

!     2nd order 
      zjt6a = zjt6a + (p(2) * za(2,n)) 

!     3rd order 
      zjt6a = zjt6a + (p(3) * za(3,n)) 

!     4th order 
      zjt6a = zjt6a + (p(4) * za(4,n) + p(5) * za(5,n)) 

!     5th order
      zjt6a = zjt6a + (p(6) * za(6,n) + p(7) * za(7,n)) 

!     6th order 
      zjt6a = zjt6a + (p(8) * za(8,n) + p(9) * za(9,n)) 

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z Jahn-Teller matrix elements up_ to 6th order
c for eb
c uses not more then 9 Parameters
c 1.: 1
c 2.: 1
c 3.: 1
c 4.: 2
c 5.: 2
c 6.: 2
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezjt6b(zjt6a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running index

      integer j               !Number of given Parameters
      integer lnx             !Maximum number of Parameters
      parameter (lnx=9)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector
 

      integer iref, ii, jj !used for ref. Hamiltonien

      double precision zjt6a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo
      
c     this takes care of split reference plus correction scheme
      ii=1
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo

!     1st order
      zjt6a = p(1) * zb(1,n)

!     2nd order 
      zjt6a = zjt6a + (p(2) * zb(2,n)) 

!     3rd order 
      zjt6a = zjt6a + (p(3) * zb(3,n)) 

!     4th order 
      zjt6a = zjt6a + (p(4) * zb(4,n) + p(5) * zb(5,n))  

!     5th order
      zjt6a = zjt6a + (p(6) * zb(6,n) + p(7) * zb(7,n))  

!     6th order 
      zjt6a = zjt6a + (p(8) * zb(8,n) + p(9) * zb(9,n))  

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z Jahn-Teller matrix elements up_ to 6th order
c for ec
c uses not more then 9 Parameters
c 1.: 1
c 2.: 1
c 3.: 1
c 4.: 2
c 5.: 2
c 6.: 2
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezjt6c(zjt6a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running index

      integer j               !Number of given Parameters
      integer lnx             !Maximum number of Parameters
      parameter (lnx=9)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector
 

      integer iref, ii, jj !used for ref. Hamiltonien

      double precision zjt6a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo
      
c     this takes care of split reference plus correction scheme
      ii=1
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo

!     1st order
      zjt6a = p(1) * zc(1,n)

!     2nd order 
      zjt6a = zjt6a + (p(2) * zc(2,n))  

!     3rd order 
      zjt6a = zjt6a + (p(3) * zc(3,n))  

!     4th order 
      zjt6a = zjt6a + (p(4) * zc(4,n) + p(5) * zc(5,n))  

!     5th order
      zjt6a = zjt6a + (p(6) * zc(6,n) + p(7) * zc(7,n))  

!     6th order 
      zjt6a = zjt6a + (p(8) * zc(8,n) + p(9) * zc(9,n))  

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Generates the diagonalelemen of the ea-eb-mode-mode terms
c up_ to sixth order
c It uses a maxiumum number of 28 parameters
c 2.:  1
c 3.:  2
c 4.:  4
c 5.:  8
c 6.: 13
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_eeab(v_eeab, n, par_, j, iref)
      implicit none

      integer i !running indices

      integer n !point
      
      integer j               !given number of parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=28)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii,jj,iref !used for Refference Hamiltonien

      double precision v_eeab !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2 !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo

!     2nd order
      v_eeab = veeab(1,n) * p(1)  

!     3rd order
      v_eeab = v_eeab + (veeab(2,n) * p(2) + veeab(3,n) * p(3))  

!     4th order
      v_eeab = v_eeab + (veeab(4,n) * p(4) + veeab(5,n) * p(5) 
     $     + veeab(6,n) * p(6) +  veeab(7,n) * p(7)) 

!     5th order
      v_eeab = v_eeab + p(8)*veeab(8,n) + p(9)*veeab(9,n)
     &     + p(10)*veeab(10,n) + p(11)*veeab(11,n) + p(12)*veeab(12,n)
     &     + p(13)*veeab(13,n) + p(14)*veeab(14,n) + p(15)*veeab(15,n)

!     6th order
      v_eeab = v_eeab + p(16)*veeab(16,n) + p(17)*veeab(17,n) 
     &     + p(18)*veeab(18,n) + p(19)*veeab(19,n) + p(20)*veeab(20,n)
     &     + p(21)*veeab(21,n) + p(22)*veeab(22,n) + p(23)*veeab(23,n)
     &     + p(24)*veeab(24,n) + p(25)*veeab(25,n) + p(26)*veeab(26,n)
     &     + p(27)*veeab(27,n) + p(28)*veeab(28,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Generates the diagonalelemen of the ea-ec-mode-mode terms
c up_ to sixth order
c It uses a maxiumum number of 28 parameters
c 2.:  1
c 3.:  2
c 4.:  4
c 5.:  8
c 6.: 13
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_eeac(v_eeab, n, par_, j, iref)
      implicit none

      integer i !running indices

      integer n !point
      
      integer j               !given number of parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=28)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii,jj,iref !used for Refference Hamiltonien

      double precision v_eeab !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo

!     2nd order
      v_eeab = veeac(1,n) * p(1) 

!     3rd order
      v_eeab = v_eeab + (veeac(2,n) * p(2) + veeac(3,n) * p(3)) 

!     4th order
      v_eeab = v_eeab + (veeac(4,n) * p(4) + veeac(5,n) * p(5) 
     $     + veeac(6,n) * p(6) +  veeac(7,n) * p(7)) 

!     5th order
      v_eeab = v_eeab + p(8)*veeac(8,n) + p(9)*veeac(9,n)
     &     + p(10)*veeac(10,n) + p(11)*veeac(11,n) + p(12)*veeac(12,n)
     &     + p(13)*veeac(13,n) + p(14)*veeac(14,n) + p(15)*veeac(15,n)

!     6th order
      v_eeab = v_eeab + p(16)*veeac(16,n) + p(17)*veeac(17,n) 
     &     + p(18)*veeac(18,n) + p(19)*veeac(19,n) + p(20)*veeac(20,n)
     &     + p(21)*veeac(21,n) + p(22)*veeac(22,n) + p(23)*veeac(23,n)
     &     + p(24)*veeac(24,n) + p(25)*veeac(25,n) + p(26)*veeac(26,n)
     &     + p(27)*veeac(27,n) + p(28)*veeac(28,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Generates the diagonalelemen of the eb-ec-mode-mode terms
c up_ to sixth order
c It uses a maxiumum number of 28 parameters
c 2.:  1
c 3.:  2
c 4.:  4
c 5.:  8
c 6.: 13
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_eebc(v_eeab, n, par_, j, iref)
      implicit none

      integer i !running indices

      integer n !point
      
      integer j               !given number of parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=28)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii,jj,iref !used for Refference Hamiltonien

      double precision v_eeab !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj 
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order
      v_eeab = veebc(1,n) * p(1) 

!     3rd order
      v_eeab = v_eeab + (veebc(2,n) * p(2) + veebc(3,n) * p(3)) 

!     4th order
      v_eeab = v_eeab + (veebc(4,n) * p(4) + veebc(5,n) * p(5) 
     $     + veebc(6,n) * p(6) +  veebc(7,n) * p(7)) 

!     5th order
      v_eeab = v_eeab + p(8)*veebc(8,n) + p(9)*veebc(9,n)
     &     + p(10)*veebc(10,n) + p(11)*veebc(11,n) + p(12)*veebc(12,n)
     &     + p(13)*veebc(13,n) + p(14)*veebc(14,n) + p(15)*veebc(15,n)

!     6th order
      v_eeab = v_eeab + p(16)*veebc(16,n) + p(17)*veebc(17,n)
     &     + p(18)*veebc(18,n) + p(19)*veebc(19,n) + p(20)*veebc(20,n)
     &     + p(21)*veebc(21,n) + p(22)*veebc(22,n) + p(23)*veebc(23,n)
     &     + p(24)*veebc(24,n) + p(25)*veebc(25,n) + p(26)*veebc(26,n)
     &     + p(27)*veebc(27,n) + p(28)*veebc(28,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Generates the mode mode coupling of the W element for ea and eb
c uses not more than 49 parameters
c 2.:  1
c 3.:  4
c 4.:  9
c 5.: 12
c 6.: 23
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_eeab(wcoup_ee, n,par_, j, iref)
      implicit none

      integer i !running indices

      integer n !point
      
      integer j               !given number of parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=49)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii,jj,iref !used for Refference Hamiltonien

      double precision wcoup_ee

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo
c     this takes care of split reference plus correction scheme
      ii=1               
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     2nd order      
      wcoup_ee = weeab(1,n) * p(1)

!     3rd order
      wcoup_ee = wcoup_ee + (weeab(2,n) * p(2) + weeab(3,n) * p(3)
     & + weeab(4,n) * p(4) + weeab(5,n) * p(5)) !/ fact(3)

!     4th order
      wcoup_ee = wcoup_ee                                         
     & +(weeab(6,n) * p( 6) + weeab(7,n) * p( 7) + weeab(8,n) * p( 8)   
     & + weeab(9,n) * p( 9) + weeab(10,n) * p(10) + weeab(11,n) * p(11) 
     & + weeab(12,n) * p(12) + weeab(13,n) * p(13) + weeab(14,n) * p(14)
     & ) !/ fact(4)                      

!     5th order
      wcoup_ee = wcoup_ee
     & +(weeab(15,n) * p(15) + weeab(16,n) * p(16) + weeab(17,n) * p(17)
     & + weeab(18,n) * p(18) + weeab(19,n) * p(19) + weeab(20,n) * p(20)
     & + weeab(21,n) * p(21) + weeab(22,n) * p(22) + weeab(23,n) * p(23)
     & + weeab(24,n)* p(24) + weeab(25,n)* p(25) + weeab(26,n)* p(26)   
     & ) !/ fact(5)                      

!     6th order
      wcoup_ee = wcoup_ee                                         
     & +(weeab(27,n) * p(27) + weeab(28,n) * p(28) + weeab(29,n) * p(29)
     & + weeab(30,n) * p(30) + weeab(31,n) * p(31) + weeab(32,n) * p(32)
     & + weeab(33,n) * p(33) + weeab(34,n) * p(34) + weeab(35,n) * p(35)
     & + weeab(36,n)* p(36) + weeab(37,n)* p(37) + weeab(38,n)* p(38)
     & + weeab(39,n)* p(39) + weeab(40,n)* p(40) + weeab(41,n)* p(41)
     & + weeab(42,n)* p(42) + weeab(43,n)* p(43) + weeab(44,n)* p(44)
     & + weeab(45,n)* p(45) + weeab(46,n)* p(46) + weeab(47,n)* p(47)
     & + weeab(48,n)* p(48) + weeab(49,n)* p(49)                   
     & ) !/ fact(6)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Generates the mode mode coupling of the W element for ea and ec
c uses not more than 49 parameters
c 2.:  1
c 3.:  4
c 4.:  9
c 5.: 12
c 6.: 23
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_eeac(wcoup_ee, n,par_, j, iref)
      implicit none

      integer i !running indices

      integer n !point
      
      integer j               !given number of parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=49)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii,jj,iref !used for Refference Hamiltonien

      double precision wcoup_ee

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1               
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order      
      wcoup_ee = weeac(1,n) * p(1)

!     3rd order
      wcoup_ee = wcoup_ee + (weeac(2,n) * p(2) + weeac(3,n) * p(3)
     & + weeac(4,n) * p(4) + weeac(5,n) * p(5)) !/ fact(3)

!     4th order
      wcoup_ee = wcoup_ee                                         
     & +(weeac(6,n) * p( 6) + weeac(7,n) * p( 7) + weeac(8,n) * p( 8)   
     & + weeac(9,n) * p( 9) + weeac(10,n) * p(10) + weeac(11,n) * p(11) 
     & + weeac(12,n) * p(12) + weeac(13,n) * p(13) + weeac(14,n) * p(14)
     & ) !/ fact(4)                      

!     5th order
      wcoup_ee = wcoup_ee
     & +(weeac(15,n) * p(15) + weeac(16,n) * p(16) + weeac(17,n) * p(17)
     & + weeac(18,n) * p(18) + weeac(19,n) * p(19) + weeac(20,n) * p(20)
     & + weeac(21,n) * p(21) + weeac(22,n) * p(22) + weeac(23,n) * p(23)
     & + weeac(24,n)* p(24) + weeac(25,n)* p(25) + weeac(26,n)* p(26)   
     & ) !/ fact(5)                      

!     6th order
      wcoup_ee = wcoup_ee                                         
     & +(weeac(27,n) * p(27) + weeac(28,n) * p(28) + weeac(29,n) * p(29)
     & + weeac(30,n) * p(30) + weeac(31,n) * p(31) + weeac(32,n) * p(32)
     & + weeac(33,n) * p(33) + weeac(34,n) * p(34) + weeac(35,n) * p(35)
     & + weeac(36,n)* p(36) + weeac(37,n)* p(37) + weeac(38,n)* p(38)
     & + weeac(39,n)* p(39) + weeac(40,n)* p(40) + weeac(41,n)* p(41)
     & + weeac(42,n)* p(42) + weeac(43,n)* p(43) + weeac(44,n)* p(44)
     & + weeac(45,n)* p(45) + weeac(46,n)* p(46) + weeac(47,n)* p(47)
     & + weeac(48,n)* p(48) + weeac(49,n)* p(49)                   
     & ) !/ fact(6)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Generates the mode mode coupling of the W element for eb and ec
c uses not more than 49 parameters
c 2.:  1
c 3.:  4
c 4.:  9
c 5.: 12
c 6.: 23
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_eebc(wcoup_ee, n,par_, j, iref)
      implicit none

      integer i !running indices

      integer n !point
      
      integer j               !given number of parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=49)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii,jj,iref !used for Refference Hamiltonien

      double precision wcoup_ee

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1               
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order      
      wcoup_ee = weebc(1,n) * p(1)

!     3rd order
      wcoup_ee = wcoup_ee + (weebc(2,n) * p(2) + weebc(3,n) * p(3)
     & + weebc(4,n) * p(4) + weebc(5,n) * p(5)) !/ fact(3)

!     4th order
      wcoup_ee = wcoup_ee                                         
     & +(weebc(6,n) * p( 6) + weebc(7,n) * p( 7) + weebc(8,n) * p( 8)   
     & + weebc(9,n) * p( 9) + weebc(10,n) * p(10) + weebc(11,n) * p(11) 
     & + weebc(12,n) * p(12) + weebc(13,n) * p(13) + weebc(14,n) * p(14)
     & ) !/ fact(4)                      

!     5th order
      wcoup_ee = wcoup_ee
     & +(weebc(15,n) * p(15) + weebc(16,n) * p(16) + weebc(17,n) * p(17)
     & + weebc(18,n) * p(18) + weebc(19,n) * p(19) + weebc(20,n) * p(20)
     & + weebc(21,n) * p(21) + weebc(22,n) * p(22) + weebc(23,n) * p(23)
     & + weebc(24,n)* p(24) + weebc(25,n)* p(25) + weebc(26,n)* p(26)   
     & ) !/ fact(5)                      

!     6th order
      wcoup_ee = wcoup_ee                                         
     & +(weebc(27,n) * p(27) + weebc(28,n) * p(28) + weebc(29,n) * p(29)
     & + weebc(30,n) * p(30) + weebc(31,n) * p(31) + weebc(32,n) * p(32)
     & + weebc(33,n) * p(33) + weebc(34,n) * p(34) + weebc(35,n) * p(35)
     & + weebc(36,n)* p(36) + weebc(37,n)* p(37) + weebc(38,n)* p(38)
     & + weebc(39,n)* p(39) + weebc(40,n)* p(40) + weebc(41,n)* p(41)
     & + weebc(42,n)* p(42) + weebc(43,n)* p(43) + weebc(44,n)* p(44)
     & + weebc(45,n)* p(45) + weebc(46,n)* p(46) + weebc(47,n)* p(47)
     & + weebc(48,n)* p(48) + weebc(49,n)* p(49)                   
     & ) !/ fact(6)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Generates the mode mode coupling of the Z element for ea and eb
c uses not more than 49 parameters
c 2.:  1
c 3.:  4
c 4.:  9
c 5.: 12
c 6.: 23
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_eeab(zcoup_ee, n,par_, j, iref)
      implicit none

      integer i !running indices

      integer n !point
      
      integer j               !given number of parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=49)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii,jj,iref !used for Refference Hamiltonien

      double precision zcoup_ee

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1 
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order      
      zcoup_ee = zeeab(1,n) * p(1)

!     3rd order
      zcoup_ee = zcoup_ee + (zeeab(2,n) * p(2) + zeeab(3,n) * p(3)
     & + zeeab(4,n) * p(4) + zeeab(5,n) * p(5)) !/ fact(3)

!     4th order
      zcoup_ee = zcoup_ee                                         
     & +(zeeab(6,n) * p( 6) + zeeab(7,n) * p( 7) + zeeab(8,n) * p( 8)   
     & + zeeab(9,n) * p( 9) + zeeab(10,n) * p(10) + zeeab(11,n) * p(11) 
     & + zeeab(12,n) * p(12) + zeeab(13,n) * p(13) + zeeab(14,n) * p(14)
     & ) !/ fact(4)                      

!     5th order
      zcoup_ee = zcoup_ee
     & +(zeeab(15,n) * p(15) + zeeab(16,n) * p(16) + zeeab(17,n) * p(17)
     & + zeeab(18,n) * p(18) + zeeab(19,n) * p(19) + zeeab(20,n) * p(20)
     & + zeeab(21,n) * p(21) + zeeab(22,n) * p(22) + zeeab(23,n) * p(23)
     & + zeeab(24,n)* p(24) + zeeab(25,n)* p(25) + zeeab(26,n)* p(26)   
     & ) !/ fact(5)                      

!     6th order
      zcoup_ee = zcoup_ee
     & +(zeeab(27,n) * p(27) + zeeab(28,n) * p(28) + zeeab(29,n) * p(29)
     & + zeeab(30,n) * p(30) + zeeab(31,n) * p(31) + zeeab(32,n) * p(32)
     & + zeeab(33,n) * p(33) + zeeab(34,n) * p(34) + zeeab(35,n) * p(35)
     & + zeeab(36,n)* p(36) + zeeab(37,n)* p(37) + zeeab(38,n)* p(38)
     & + zeeab(39,n)* p(39) + zeeab(40,n)* p(40) + zeeab(41,n)* p(41)
     & + zeeab(42,n)* p(42) + zeeab(43,n)* p(43) + zeeab(44,n)* p(44)
     & + zeeab(45,n)* p(45) + zeeab(46,n)* p(46) + zeeab(47,n)* p(47)
     & + zeeab(48,n)* p(48) + zeeab(49,n)* p(49)
     & ) !/ fact(6)                      

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Generates the mode mode coupling of the Z element for ea and ec
c uses not more than 49 parameters
c 2.:  1
c 3.:  4
c 4.:  9
c 5.: 12
c 6.: 23
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_eeac(zcoup_ee, n,par_, j, iref)
      implicit none

      integer i !running indices

      integer n !point
      
      integer j               !given number of parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=49)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii,jj,iref !used for Refference Hamiltonien

      double precision zcoup_ee

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1 
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order      
      zcoup_ee = zeeac(1,n) * p(1)

!     3rd order
      zcoup_ee = zcoup_ee + (zeeac(2,n) * p(2) + zeeac(3,n) * p(3)
     & + zeeac(4,n) * p(4) + zeeac(5,n) * p(5)) !/ fact(3)

!     4th order
      zcoup_ee = zcoup_ee                                         
     & +(zeeac(6,n) * p( 6) + zeeac(7,n) * p( 7) + zeeac(8,n) * p( 8)   
     & + zeeac(9,n) * p( 9) + zeeac(10,n) * p(10) + zeeac(11,n) * p(11) 
     & + zeeac(12,n) * p(12) + zeeac(13,n) * p(13) + zeeac(14,n) * p(14)
     & ) !/ fact(4)                      

!     5th order
      zcoup_ee = zcoup_ee
     & +(zeeac(15,n) * p(15) + zeeac(16,n) * p(16) + zeeac(17,n) * p(17)
     & + zeeac(18,n) * p(18) + zeeac(19,n) * p(19) + zeeac(20,n) * p(20)
     & + zeeac(21,n) * p(21) + zeeac(22,n) * p(22) + zeeac(23,n) * p(23)
     & + zeeac(24,n)* p(24) + zeeac(25,n)* p(25) + zeeac(26,n)* p(26)   
     & ) !/ fact(5)                      

!     6th order
      zcoup_ee = zcoup_ee
     & +(zeeac(27,n) * p(27) + zeeac(28,n) * p(28) + zeeac(29,n) * p(29)
     & + zeeac(30,n) * p(30) + zeeac(31,n) * p(31) + zeeac(32,n) * p(32)
     & + zeeac(33,n) * p(33) + zeeac(34,n) * p(34) + zeeac(35,n) * p(35)
     & + zeeac(36,n)* p(36) + zeeac(37,n)* p(37) + zeeac(38,n)* p(38)
     & + zeeac(39,n)* p(39) + zeeac(40,n)* p(40) + zeeac(41,n)* p(41)
     & + zeeac(42,n)* p(42) + zeeac(43,n)* p(43) + zeeac(44,n)* p(44)
     & + zeeac(45,n)* p(45) + zeeac(46,n)* p(46) + zeeac(47,n)* p(47)
     & + zeeac(48,n)* p(48) + zeeac(49,n)* p(49)
     & ) !/ fact(6)                      

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c Generates the mode mode coupling of the Z element for eb and ec
c uses not more than 49 parameters
c 2.:  1
c 3.:  4
c 4.:  9
c 5.: 12
c 6.: 23
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_eebc(zcoup_ee, n,par_, j, iref)
      implicit none

      integer i !running indices

      integer n !point
      
      integer j               !given number of parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=49)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector


      integer ii,jj,iref !used for Refference Hamiltonien

      double precision zcoup_ee

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1 
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order      
      zcoup_ee = zeebc(1,n) * p(1)

!     3rd order
      zcoup_ee = zcoup_ee + (zeebc(2,n) * p(2) + zeebc(3,n) * p(3)
     & + zeebc(4,n) * p(4) + zeebc(5,n) * p(5)) !/ fact(3)

!     4th order
      zcoup_ee = zcoup_ee                                         
     & +(zeebc(6,n) * p( 6) + zeebc(7,n) * p( 7) + zeebc(8,n) * p( 8)   
     & + zeebc(9,n) * p( 9) + zeebc(10,n) * p(10) + zeebc(11,n) * p(11) 
     & + zeebc(12,n) * p(12) + zeebc(13,n) * p(13) + zeebc(14,n) * p(14)
     & ) !/ fact(4)                      

!     5th order
      zcoup_ee = zcoup_ee
     & +(zeebc(15,n) * p(15) + zeebc(16,n) * p(16) + zeebc(17,n) * p(17)
     & + zeebc(18,n) * p(18) + zeebc(19,n) * p(19) + zeebc(20,n) * p(20)
     & + zeebc(21,n) * p(21) + zeebc(22,n) * p(22) + zeebc(23,n) * p(23)
     & + zeebc(24,n)* p(24) + zeebc(25,n)* p(25) + zeebc(26,n)* p(26)   
     & ) !/ fact(5)                      

!     6th order
      zcoup_ee = zcoup_ee
     & +(zeebc(27,n) * p(27) + zeebc(28,n) * p(28) + zeebc(29,n) * p(29)
     & + zeebc(30,n) * p(30) + zeebc(31,n) * p(31) + zeebc(32,n) * p(32)
     & + zeebc(33,n) * p(33) + zeebc(34,n) * p(34) + zeebc(35,n) * p(35)
     & + zeebc(36,n)* p(36) + zeebc(37,n)* p(37) + zeebc(38,n)* p(38)
     & + zeebc(39,n)* p(39) + zeebc(40,n)* p(40) + zeebc(41,n)* p(41)
     & + zeebc(42,n)* p(42) + zeebc(43,n)* p(43) + zeebc(44,n)* p(44)
     & + zeebc(45,n)* p(45) + zeebc(46,n)* p(46) + zeebc(47,n)* p(47)
     & + zeebc(48,n)* p(48) + zeebc(49,n)* p(49)
     & ) !/ fact(6)                      

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling betweeabn ea
c and a modes up_ to sixth order
c uses not more than 10 parameters
c 3.: 1
c 4.: 2
c 5.: 3
c 6.: 4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_eaa(v_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=10)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision v_e1a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo
!.....
!     3rd order
      v_e1a =  (p(1)*pa(1,n)*va(1,n))
 
!     4th order
      v_e1a = v_e1a+(p(2)*pa(1,n)*va(2,n)+p(3)*pa(2,n)*va(1,n))

!     5th order
      v_e1a = v_e1a + (p(4)*pa(1,n)*va(3,n)+p(5)*pa(2,n)*va(2,n)
     &        + p(6)*pa(3,n)*va(1,n))

!     6th order
      v_e1a= v_e1a + p(7)*pa(1,n)*va(4,n) + p(8)*pa(2,n)*va(3,n)
     &        + p(9)*pa(3,n)*va(2,n) + p(10)*pa(4,n)*va(1,n)
      
      end      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling betweeabn eb
c and a modes up_ to sixth order
c uses not more than 10 parameters
c 3.: 1
c 4.: 2
c 5.: 3
c 6.: 4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_eba(v_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=10)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision v_e1a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo
!.....
!     3rd order
      v_e1a =  (p(1)*pa(1,n)*vb(1,n))
 
!     4th order
      v_e1a = v_e1a+(p(2)*pa(1,n)*vb(2,n)+p(3)*pa(2,n)*vb(1,n))

!     5th order
      v_e1a = v_e1a + (p(4)*pa(1,n)*vb(3,n)+p(5)*pa(2,n)*vb(2,n)
     &        + p(6)*pa(3,n)*vb(1,n))

!     6th order
      v_e1a= v_e1a + p(7)*pa(1,n)*vb(4,n) + p(8)*pa(2,n)*vb(3,n)
     &        + p(9)*pa(3,n)*vb(2,n) + p(10)*pa(4,n)*vb(1,n)
      
      end      

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling betweeabn ec
c and ac modes up_ to sixth order
c uses not more than 10 parameters
c 3.: 1
c 4.: 2
c 5.: 3
c 6.: 4 
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_eca(v_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=10)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision v_e1a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

      ! this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo
!.....
!     3rd order
      v_e1a =  (p(1)*pa(1,n)*vc(1,n))
 
!     4th order
      v_e1a = v_e1a+(p(2)*pa(1,n)*vc(2,n)+p(3)*pa(2,n)*vc(1,n))

!     5th order
      v_e1a = v_e1a + (p(4)*pa(1,n)*vc(3,n)+p(5)*pa(2,n)*vc(2,n)
     &        + p(6)*pa(3,n)*vc(1,n))

!     6th order
      v_e1a= v_e1a + p(7)*pa(1,n)*vc(4,n) + p(8)*pa(2,n)*vc(3,n)
     &        + p(9)*pa(3,n)*vc(2,n) + p(10)*pa(4,n)*vc(1,n)
      
      end      

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling betweeabn
c ea and b modes up_ to sixth order
c uses not more than 10 parameters
c 3.: 1
c 4.: 2
c 5.: 3
c 6.: 4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_eab(v_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=10)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision v_e1a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo
!.....
!     3rd order
      v_e1a =  (p(1)*pb(1,n)*va(1,n))

!     4th order
      v_e1a = v_e1a+(p(2)*pb(1,n)*va(2,n)+p(3)*pb(2,n)*va(1,n))

!     5th order
      v_e1a = v_e1a + (p(4)*pb(1,n)*va(3,n)+p(5)*pb(2,n)*va(2,n)
     &        + p(6)*pb(3,n)*va(1,n))

!     6th order
      v_e1a= v_e1a + p(7)*pb(1,n)*va(4,n) + p(8)*pb(2,n)*va(3,n)
     &        + p(9)*pb(3,n)*va(2,n) + p(10)*pb(4,n)*va(1,n)
      
      end      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling betweeabn e
c b and b modes up_ to sixth order
c uses not more than 10 parameters
c 3.: 1
c 4.: 2
c 5.: 3
c 6.: 4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_ebb(v_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=10)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision v_e1a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo
!.....
!     3rd order
      v_e1a =  (p(1)*pb(1,n)*vb(1,n))
 
!     4th order
      v_e1a = v_e1a+(p(2)*pb(1,n)*vb(2,n)+p(3)*pb(2,n)*vb(1,n))

!     5th order
      v_e1a = v_e1a + (p(4)*pb(1,n)*vb(3,n)+p(5)*pb(2,n)*vb(2,n)
     &        + p(6)*pb(3,n)*vb(1,n))

!     6th order
      v_e1a= v_e1a + p(7)*pb(1,n)*vb(4,n) + p(8)*pb(2,n)*vb(3,n)
     &        + p(9)*pb(3,n)*vb(2,n) + p(10)*pb(4,n)*vb(1,n)
      
      end      

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling betweeabn ec
c and b modes up_ to sixth order
c uses not more than 10 parameters
c 3.: 1
c 4.: 2
c 5.: 3
c 6.: 4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_ecb(v_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=10)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision v_e1a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo
!.....
!     3rd order
      v_e1a =  (p(1)*pb(1,n)*vc(1,n))
 
!     4th order
      v_e1a = v_e1a+(p(2)*pb(1,n)*vc(2,n)+p(3)*pb(2,n)*vc(1,n))

!     5th order
      v_e1a = v_e1a + (p(4)*pb(1,n)*vc(3,n)+p(5)*pb(2,n)*vc(2,n)
     &        + p(6)*pb(3,n)*vc(1,n))

!     6th order
      v_e1a= v_e1a + p(7)*pb(1,n)*vc(4,n) + p(8)*pb(2,n)*vc(3,n)
     &        + p(9)*pb(3,n)*vc(2,n) + p(10)*pb(4,n)*vc(1,n)
      
      end      

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling betweeabn ea
c and modes up_ to sixth order
c uses not more than 10 parameters
c 3.: 1
c 4.: 2
c 5.: 3
c 6.: 4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_eac(v_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=10)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision v_e1a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo
!.....
!     3rd order
      v_e1a =  (p(1)*pc(1,n)*va(1,n))
 
!     4th order
      v_e1a = v_e1a+(p(2)*pc(1,n)*va(2,n)+p(3)*pc(2,n)*va(1,n))

!     5th order
      v_e1a = v_e1a + (p(4)*pc(1,n)*va(3,n)+p(5)*pc(2,n)*va(2,n)
     &        + p(6)*pc(3,n)*va(1,n))

!     6th order
      v_e1a= v_e1a + p(7)*pc(1,n)*va(4,n) + p(8)*pc(2,n)*va(3,n)
     &        + p(9)*pc(3,n)*va(2,n) + p(10)*pc(4,n)*va(1,n)
      
      end      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling betweeabn eb
c and modes up_ to sixth order
c uses not more than 10 parameters
c 3.: 1
c 4.: 2
c 5.: 3
c 6.: 4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_ebc(v_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=10)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision v_e1a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

C     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo
!.....
!     3rd order
      v_e1a =  (p(1)*pc(1,n)*vb(1,n))
 
!     4th order
      v_e1a = v_e1a+(p(2)*pc(1,n)*vb(2,n)+p(3)*pc(2,n)*vb(1,n))

!     5th order
      v_e1a = v_e1a + (p(4)*pc(1,n)*vb(3,n)+p(5)*pc(2,n)*vb(2,n)
     &        + p(6)*pc(3,n)*vb(1,n))

!     6th order
      v_e1a= v_e1a + p(7)*pc(1,n)*vb(4,n) + p(8)*pc(2,n)*vb(3,n)
     &        + p(9)*pc(3,n)*vb(2,n) + p(10)*pc(4,n)*vb(1,n)
      
      end      

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling betweeabn ec
C and modes up_ to sixth order
c uses not more than 10 parameters
c 3.: 1
c 4.: 2
c 5.: 3
c 6.: 4
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_ecc(v_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=10)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision v_e1a !Matrix element
      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

      ! this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo
!.....
!     3rd order
      v_e1a =  (p(1)*pc(1,n)*vc(1,n))
 
!     4th order
      v_e1a = v_e1a+(p(2)*pc(1,n)*vc(2,n)+p(3)*pc(2,n)*vc(1,n))

!     5th order
      v_e1a = v_e1a + (p(4)*pc(1,n)*vc(3,n)+p(5)*pc(2,n)*vc(2,n)
     &        + p(6)*pc(3,n)*vc(1,n))

!     6th order
      v_e1a= v_e1a + p(7)*pc(1,n)*vc(4,n) + p(8)*pc(2,n)*vc(3,n)
     &        + p(9)*pc(3,n)*vc(2,n) + p(10)*pc(4,n)*vc(1,n)
      
      end      

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between ea and a
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_eaa(wcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision wcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2   !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      wcoup_e1a = p(1)*pa(1,n)*wa(1,n)

!     3rd order 
      wcoup_e1a = wcoup_e1a+(p(2)*pa(1,n)*wa(2,n)+p(3)*pa(2,n)*wa(1,n))

!     4th order 
      wcoup_e1a = wcoup_e1a+(p(4)*pa(1,n)*wa(3,n)+p(5)*pa(2,n)*wa(2,n)
     & + p(6)*pa(3,n)*wa(1,n))

!     5th order 
      wcoup_e1a = wcoup_e1a + (p(7)*pa(4,n)*wa(1,n)+p(8)*pa(3,n)*wa(2,n)
     & +p(9)*pa(2,n)*wa(3,n)+p(10)*pa(1,n)*wa(4,n)
     & +p(11)*pa(1,n)*wa(5,n))

!     6th order 
      wcoup_e1a = wcoup_e1a + 
     &    p(12)*pa(5,n)*wa(1,n)+p(13)*pa(4,n)*wa(2,n)
     &    + p(14)*pa(3,n)*wa(3,n)
     &    + p(15)*pa(2,n)*wa(4,n)+p(16)*pa(2,n)*wa(5,n)
     &    + pa(1,n)*(p(17)*wa(6,n)+p(18)*wa(7,n))

      end      
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between eb and a
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_eba(wcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision wcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

      ! this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2   !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      wcoup_e1a = p(1)*pa(1,n)*wb(1,n)

!     3rd order 
      wcoup_e1a = wcoup_e1a+(p(2)*pa(1,n)*wb(2,n)+p(3)*pa(2,n)*wb(1,n))

!     4th order 
      wcoup_e1a = wcoup_e1a+(p(4)*pa(1,n)*wb(3,n)+p(5)*pa(2,n)*wb(2,n)
     & + p(6)*pa(3,n)*wb(1,n))

!     5th order 
      wcoup_e1a = wcoup_e1a + (p(7)*pa(4,n)*wb(1,n)+p(8)*pa(3,n)*wb(2,n)
     & +p(9)*pa(2,n)*wb(3,n)+p(10)*pa(1,n)*wb(4,n)
     & +p(11)*pa(1,n)*wb(5,n))

!     6th order 
      wcoup_e1a = wcoup_e1a + 
     &    p(12)*pa(5,n)*wb(1,n)+p(13)*pa(4,n)*wb(2,n)
     &    + p(14)*pa(3,n)*wb(3,n)
     &    + p(15)*pa(2,n)*wb(4,n)+p(16)*pa(2,n)*wb(5,n)
     &    + pa(1,n)*(p(17)*wb(6,n)+p(18)*wb(7,n))

      end      
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between ec and a
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_eca(wcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision wcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

      ! this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2   !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      wcoup_e1a = p(1)*pa(1,n)*wc(1,n)

!     3rd order 
      wcoup_e1a = wcoup_e1a+(p(2)*pa(1,n)*wc(2,n)+p(3)*pa(2,n)*wc(1,n))

!     4th order 
      wcoup_e1a = wcoup_e1a+(p(4)*pa(1,n)*wc(3,n)+p(5)*pa(2,n)*wc(2,n)
     & + p(6)*pa(3,n)*wc(1,n))

!     5th order 
      wcoup_e1a = wcoup_e1a + (p(7)*pa(4,n)*wc(1,n)+p(8)*pa(3,n)*wc(2,n)
     & +p(9)*pa(2,n)*wc(3,n)+p(10)*pa(1,n)*wc(4,n)
     & +p(11)*pa(1,n)*wc(5,n))

!     6th order 
      wcoup_e1a = wcoup_e1a + 
     &    p(12)*pa(5,n)*wc(1,n)+p(13)*pa(4,n)*wc(2,n)
     &    + p(14)*pa(3,n)*wc(3,n)
     &    + p(15)*pa(2,n)*wc(4,n)+p(16)*pa(2,n)*wc(5,n)
     &    + pa(1,n)*(p(17)*wc(6,n)+p(18)*wc(7,n))

      end
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between ea and b
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_eab(wcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision wcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2   !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      wcoup_e1a = p(1)*pb(1,n)*wa(1,n)

!     3rd order 
      wcoup_e1a = wcoup_e1a+(p(2)*pb(1,n)*wa(2,n)+p(3)*pb(2,n)*wa(1,n))        

!     4th order 
      wcoup_e1a = wcoup_e1a+(p(4)*pb(1,n)*wa(3,n)+p(5)*pb(2,n)*wa(2,n)
     & + p(6)*pb(3,n)*wa(1,n))

!     5th order 
      wcoup_e1a = wcoup_e1a + (p(7)*pb(4,n)*wa(1,n)+p(8)*pb(3,n)*wa(2,n)
     & +p(9)*pb(2,n)*wa(3,n)+p(10)*pb(1,n)*wa(4,n)
     & +p(11)*pb(1,n)*wa(5,n))

!     6th order 
      wcoup_e1a = wcoup_e1a + 
     &    p(12)*pb(5,n)*wa(1,n)+p(13)*pb(4,n)*wa(2,n)
     &    + p(14)*pb(3,n)*wa(3,n)
     &    + p(15)*pb(2,n)*wa(4,n)+p(16)*pb(2,n)*wa(5,n)
     &    + pb(1,n)*(p(17)*wa(6,n)+p(18)*wa(7,n))

      end      
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between eb and b
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_ebb(wcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision wcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2   !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      wcoup_e1a = p(1)*pb(1,n)*wb(1,n)

!     3rd order 
      wcoup_e1a = wcoup_e1a+(p(2)*pb(1,n)*wb(2,n)+p(3)*pb(2,n)*wb(1,n))

!     4th order 
      wcoup_e1a = wcoup_e1a+(p(4)*pb(1,n)*wb(3,n)+p(5)*pb(2,n)*wb(2,n)
     & + p(6)*pb(3,n)*wb(1,n))

!     5th order 
      wcoup_e1a = wcoup_e1a + (p(7)*pb(4,n)*wb(1,n)+p(8)*pb(3,n)*wb(2,n)
     & +p(9)*pb(2,n)*wb(3,n)+p(10)*pb(1,n)*wb(4,n)
     & +p(11)*pb(1,n)*wb(5,n))

!     6th order 
      wcoup_e1a = wcoup_e1a + 
     &    p(12)*pb(5,n)*wb(1,n)+p(13)*pb(4,n)*wb(2,n)
     &    + p(14)*pb(3,n)*wb(3,n)
     &    + p(15)*pb(2,n)*wb(4,n)+p(16)*pb(2,n)*wb(5,n)
     &    + pb(1,n)*(p(17)*wb(6,n)+p(18)*wb(7,n))

      end      
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between ec and b
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_ecb(wcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision wcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2   !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      wcoup_e1a = p(1)*pb(1,n)*wc(1,n)

!     3rd order 
      wcoup_e1a = wcoup_e1a+(p(2)*pb(1,n)*wc(2,n)+p(3)*pb(2,n)*wc(1,n))

!     4th order 
      wcoup_e1a = wcoup_e1a+(p(4)*pb(1,n)*wc(3,n)+p(5)*pb(2,n)*wc(2,n)
     & + p(6)*pb(3,n)*wc(1,n))

!     5th order 
      wcoup_e1a = wcoup_e1a + (p(7)*pb(4,n)*wc(1,n)+p(8)*pb(3,n)*wc(2,n)
     & +p(9)*pb(2,n)*wc(3,n)+p(10)*pb(1,n)*wc(4,n)
     & +p(11)*pb(1,n)*wc(5,n))

!     6th order 
      wcoup_e1a = wcoup_e1a + 
     &    p(12)*pb(5,n)*wc(1,n)+p(13)*pb(4,n)*wc(2,n)
     &    + p(14)*pb(3,n)*wc(3,n)
     &    + p(15)*pb(2,n)*wc(4,n)+p(16)*pb(2,n)*wc(5,n)
     &    + pb(1,n)*(p(17)*wc(6,n)+p(18)*wc(7,n))

      end
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between ea and c
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_eac(wcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision wcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2   !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      wcoup_e1a = p(1)*pc(1,n)*wa(1,n)

!     3rd order 
      wcoup_e1a = wcoup_e1a+(p(2)*pc(1,n)*wa(2,n)+p(3)*pc(2,n)*wa(1,n))

!     4th order 
      wcoup_e1a = wcoup_e1a+(p(4)*pc(1,n)*wa(3,n)+p(5)*pc(2,n)*wa(2,n)
     & + p(6)*pc(3,n)*wa(1,n))

!     5th order 
      wcoup_e1a = wcoup_e1a + (p(7)*pc(4,n)*wa(1,n)+p(8)*pc(3,n)*wa(2,n)
     & +p(9)*pc(2,n)*wa(3,n)+p(10)*pc(1,n)*wa(4,n)
     & +p(11)*pc(1,n)*wa(5,n))

!     6th order 
      wcoup_e1a = wcoup_e1a + 
     &    p(12)*pc(5,n)*wa(1,n)+p(13)*pc(4,n)*wa(2,n)
     &    + p(14)*pc(3,n)*wa(3,n)
     &    + p(15)*pc(2,n)*wa(4,n)+p(16)*pc(2,n)*wa(5,n)
     &    + pc(1,n)*(p(17)*wa(6,n)+p(18)*wa(7,n))

      end      
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between eb and c
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_ebc(wcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision wcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2   !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      wcoup_e1a = p(1)*pc(1,n)*wb(1,n)

!     3rd order 
      wcoup_e1a = wcoup_e1a+(p(2)*pc(1,n)*wb(2,n)+p(3)*pc(2,n)*wb(1,n))

!     4th order 
      wcoup_e1a = wcoup_e1a+(p(4)*pc(1,n)*wb(3,n)+p(5)*pc(2,n)*wb(2,n)
     & + p(6)*pc(3,n)*wb(1,n))

!     5th order 
      wcoup_e1a = wcoup_e1a + (p(7)*pc(4,n)*wb(1,n)+p(8)*pc(3,n)*wb(2,n)
     & +p(9)*pc(2,n)*wb(3,n)+p(10)*pc(1,n)*wb(4,n)
     & +p(11)*pc(1,n)*wb(5,n))

!     6th order 
      wcoup_e1a = wcoup_e1a + 
     &    p(12)*pc(5,n)*wb(1,n)+p(13)*pc(4,n)*wb(2,n)
     &    + p(14)*pc(3,n)*wb(3,n)
     &    + p(15)*pc(2,n)*wb(4,n)+p(16)*pc(2,n)*wb(5,n)
     &    + pc(1,n)*(p(17)*wb(6,n)+p(18)*wb(7,n))

      end      
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between ec and c
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ewcoup_ecc(wcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision wcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2   !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      wcoup_e1a = p(1)*pc(1,n)*wc(1,n)

!     3rd order 
      wcoup_e1a = wcoup_e1a+(p(2)*pc(1,n)*wc(2,n)+p(3)*pc(2,n)*wc(1,n))

!     4th order 
      wcoup_e1a = wcoup_e1a+(p(4)*pc(1,n)*wc(3,n)+p(5)*pc(2,n)*wc(2,n)
     & + p(6)*pc(3,n)*wc(1,n))

!     5th order 
      wcoup_e1a = wcoup_e1a + (p(7)*pc(4,n)*wc(1,n)+p(8)*pc(3,n)*wc(2,n)
     & +p(9)*pc(2,n)*wc(3,n)+p(10)*pc(1,n)*wc(4,n)
     & +p(11)*pc(1,n)*wc(5,n))

!     6th order 
      wcoup_e1a = wcoup_e1a + 
     &    p(12)*pc(5,n)*wc(1,n)+p(13)*pc(4,n)*wc(2,n)
     &    + p(14)*pc(3,n)*wc(3,n)
     &    + p(15)*pc(2,n)*wc(4,n)+p(16)*pc(2,n)*wc(5,n)
     &    + pc(1,n)*(p(17)*wc(6,n)+p(18)*wc(7,n))

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between ea and a
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_eaa(zcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision zcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2   !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      zcoup_e1a = p(1)*pa(1,n)*za(1,n)

!     3rd order 
      zcoup_e1a=zcoup_e1a+ (p(2)*pa(1,n)*za(2,n) + p(3)*pa(2,n)*za(1,n))

!     4th order 
      zcoup_e1a=zcoup_e1a + (p(4)*pa(1,n)*za(3,n) + p(5)*pa(2,n)*za(2,n)
     & + p(6)*pa(3,n)*za(1,n))

!     5th order 
      zcoup_e1a=zcoup_e1a + (p(7)*pa(4,n)*za(1,n) + p(8)*pa(3,n)*za(2,n)
     & +p(9)*pa(2,n)*za(3,n)+p(10)*pa(1,n)*za(4,n)
     & +p(11)*pa(1,n)*za(5,n))

!     6th order 
      zcoup_e1a = zcoup_e1a + 
     &    p(12)*pa(5,n)*za(1,n)+p(13)*pa(4,n)*za(2,n)
     &   +p(14)*pa(3,n)*za(3,n)
     &    + p(15)*pa(2,n)*za(4,n)+p(16)*pa(2,n)*za(5,n)
     &    + pa(1,n)*(p(17)*za(6,n)+p(18)*za(7,n))


      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between eb and a
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_eba(zcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision zcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2   !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      zcoup_e1a = p(1)*pa(1,n)*zb(1,n)

!     3rd order 
      zcoup_e1a=zcoup_e1a+ (p(2)*pa(1,n)*zb(2,n) + p(3)*pa(2,n)*zb(1,n))

!     4th order 
      zcoup_e1a=zcoup_e1a + (p(4)*pa(1,n)*zb(3,n) + p(5)*pa(2,n)*zb(2,n)
     & + p(6)*pa(3,n)*zb(1,n))

!     5th order 
      zcoup_e1a=zcoup_e1a + (p(7)*pa(4,n)*zb(1,n) + p(8)*pa(3,n)*zb(2,n)
     & +p(9)*pa(2,n)*zb(3,n)+p(10)*pa(1,n)*zb(4,n)
     & +p(11)*pa(1,n)*zb(5,n))

!     6th order 
      zcoup_e1a = zcoup_e1a + 
     &    p(12)*pa(5,n)*zb(1,n)+p(13)*pa(4,n)*zb(2,n)
     &   +p(14)*pa(3,n)*zb(3,n)
     &    + p(15)*pa(2,n)*zb(4,n)+p(16)*pa(2,n)*zb(5,n)
     &    + pa(1,n)*(p(17)*zb(6,n)+p(18)*zb(7,n))


      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between ec and a
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_eca(zcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision zcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      zcoup_e1a = p(1)*pa(1,n)*zc(1,n)

!     3rd order 
      zcoup_e1a=zcoup_e1a+ (p(2)*pa(1,n)*zc(2,n) + p(3)*pa(2,n)*zc(1,n))

!     4th order 
      zcoup_e1a=zcoup_e1a + (p(4)*pa(1,n)*zc(3,n) + p(5)*pa(2,n)*zc(2,n)
     & + p(6)*pa(3,n)*zc(1,n)) 

!     5th order 
      zcoup_e1a=zcoup_e1a + (p(7)*pa(4,n)*zc(1,n) + p(8)*pa(3,n)*zc(2,n)
     & +p(9)*pa(2,n)*zc(3,n)+p(10)*pa(1,n)*zc(4,n)
     & +p(11)*pa(1,n)*zc(5,n)) 

!     6th order 
      zcoup_e1a = zcoup_e1a + 
     &    p(12)*pa(5,n)*zc(1,n)+p(13)*pa(4,n)*zc(2,n)
     &   +p(14)*pa(3,n)*zc(3,n)
     &    + p(15)*pa(2,n)*zc(4,n)+p(16)*pa(2,n)*zc(5,n)
     &    + pa(1,n)*(p(17)*zc(6,n)+p(18)*zc(7,n))


      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between ea and b
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_eab(zcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision zcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      zcoup_e1a = p(1)*pb(1,n)*za(1,n)

!     3rd order 
      zcoup_e1a=zcoup_e1a+ (p(2)*pb(1,n)*za(2,n) + p(3)*pb(2,n)*za(1,n))

!     4th order 
      zcoup_e1a=zcoup_e1a + (p(4)*pb(1,n)*za(3,n) + p(5)*pb(2,n)*za(2,n)
     & + p(6)*pb(3,n)*za(1,n))

!     5th order 
      zcoup_e1a=zcoup_e1a + (p(7)*pb(4,n)*za(1,n) + p(8)*pb(3,n)*za(2,n)
     & +p(9)*pb(2,n)*za(3,n)+p(10)*pb(1,n)*za(4,n)
     & +p(11)*pb(1,n)*za(5,n))

!     6th order 
      zcoup_e1a = zcoup_e1a + 
     &    p(12)*pb(5,n)*za(1,n)+p(13)*pb(4,n)*za(2,n)
     &   +p(14)*pb(3,n)*za(3,n)
     &    + p(15)*pb(2,n)*za(4,n)+p(16)*pb(2,n)*za(5,n)
     &    + pb(1,n)*(p(17)*za(6,n)+p(18)*za(7,n))


      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between eb and b
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_ebb(zcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision zcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      zcoup_e1a = p(1)*pb(1,n)*zb(1,n)

!     3rd order 
      zcoup_e1a=zcoup_e1a+ (p(2)*pb(1,n)*zb(2,n) + p(3)*pb(2,n)*zb(1,n))

!     4th order 
      zcoup_e1a=zcoup_e1a + (p(4)*pb(1,n)*zb(3,n) + p(5)*pb(2,n)*zb(2,n)
     & + p(6)*pb(3,n)*zb(1,n))

!     5th order 
      zcoup_e1a=zcoup_e1a + (p(7)*pb(4,n)*zb(1,n) + p(8)*pb(3,n)*zb(2,n)
     & +p(9)*pb(2,n)*zb(3,n)+p(10)*pb(1,n)*zb(4,n)
     & +p(11)*pb(1,n)*zb(5,n))

!     6th order 
      zcoup_e1a = zcoup_e1a + 
     &    p(12)*pb(5,n)*zb(1,n)+p(13)*pb(4,n)*zb(2,n)
     &   +p(14)*pb(3,n)*zb(3,n)
     &    + p(15)*pb(2,n)*zb(4,n)+p(16)*pb(2,n)*zb(5,n)
     &    + pb(1,n)*(p(17)*zb(6,n)+p(18)*zb(7,n))


      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between ec and b
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_ecb(zcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision zcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      zcoup_e1a = p(1)*pb(1,n)*zc(1,n)

!     3rd order 
      zcoup_e1a=zcoup_e1a+ (p(2)*pb(1,n)*zc(2,n) + p(3)*pb(2,n)*zc(1,n))

!     4th order 
      zcoup_e1a=zcoup_e1a + (p(4)*pb(1,n)*zc(3,n) + p(5)*pb(2,n)*zc(2,n)
     & + p(6)*pb(3,n)*zc(1,n))

!     5th order 
      zcoup_e1a=zcoup_e1a + (p(7)*pb(4,n)*zc(1,n) + p(8)*pb(3,n)*zc(2,n)
     & +p(9)*pb(2,n)*zc(3,n)+p(10)*pb(1,n)*zc(4,n)
     & +p(11)*pb(1,n)*zc(5,n))

!     6th order 
      zcoup_e1a = zcoup_e1a + 
     &    p(12)*pb(5,n)*zc(1,n)+p(13)*pb(4,n)*zc(2,n)
     &   +p(14)*pb(3,n)*zc(3,n)
     &    + p(15)*pb(2,n)*zc(4,n)+p(16)*pb(2,n)*zc(5,n)
     &    + pb(1,n)*(p(17)*zc(6,n)+p(18)*zc(7,n))


      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between ea and c
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_eac(zcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision zcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      zcoup_e1a = p(1)*pc(1,n)*za(1,n)

!     3rd order 
      zcoup_e1a=zcoup_e1a+ (p(2)*pc(1,n)*za(2,n) + p(3)*pc(2,n)*za(1,n))

!     4th order 
      zcoup_e1a=zcoup_e1a + (p(4)*pc(1,n)*za(3,n) + p(5)*pc(2,n)*za(2,n)
     & + p(6)*pc(3,n)*za(1,n))

!     5th order 
      zcoup_e1a=zcoup_e1a + (p(7)*pc(4,n)*za(1,n) + p(8)*pc(3,n)*za(2,n)
     & +p(9)*pc(2,n)*za(3,n)+p(10)*pc(1,n)*za(4,n)
     & +p(11)*pc(1,n)*za(5,n))

!     6th order 
      zcoup_e1a = zcoup_e1a + 
     &    p(12)*pc(5,n)*za(1,n)+p(13)*pc(4,n)*za(2,n)
     &   +p(14)*pc(3,n)*za(3,n)
     &    + p(15)*pc(2,n)*za(4,n)+p(16)*pc(2,n)*za(5,n)
     &    + pc(1,n)*(p(17)*za(6,n)+p(18)*za(7,n))


      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between eb and c
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_ebc(zcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision zcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      zcoup_e1a = p(1)*pc(1,n)*zb(1,n)

!     3rd order 
      zcoup_e1a=zcoup_e1a+ (p(2)*pc(1,n)*zb(2,n) + p(3)*pc(2,n)*zb(1,n))

!     4th order 
      zcoup_e1a=zcoup_e1a + (p(4)*pc(1,n)*zb(3,n) + p(5)*pc(2,n)*zb(2,n)
     & + p(6)*pc(3,n)*zb(1,n))

!     5th order 
      zcoup_e1a=zcoup_e1a + (p(7)*pc(4,n)*zb(1,n) + p(8)*pc(3,n)*zb(2,n)
     & +p(9)*pc(2,n)*zb(3,n)+p(10)*pc(1,n)*zb(4,n)
     & +p(11)*pc(1,n)*zb(5,n))

!     6th order 
      zcoup_e1a = zcoup_e1a + 
     &    p(12)*pc(5,n)*zb(1,n)+p(13)*pc(4,n)*zb(2,n)
     &   +p(14)*pc(3,n)*zb(3,n)
     &    + p(15)*pc(2,n)*zb(4,n)+p(16)*pc(2,n)*zb(5,n)
     &    + pc(1,n)*(p(17)*zb(6,n)+p(18)*zb(7,n))


      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between ec and c
c modes up_ to sixth order
c uses not more than 18 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 5
c 6.: 7
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ezcoup_ecc(zcoup_e1a, n, par_, j, iref)
      implicit none

      integer n !point

      integer i !running idices

      integer j               !number of given parameters
      integer lnx             !maximal number of parameters
      parameter (lnx=18)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !inernal parameter vector


      integer iref, ii, jj !borders for refference model

      double precision zcoup_e1a !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!.....
!     2nd order 
      zcoup_e1a = p(1)*pc(1,n)*zc(1,n)

!     3rd order 
      zcoup_e1a=zcoup_e1a+ (p(2)*pc(1,n)*zc(2,n) + p(3)*pc(2,n)*zc(1,n))

!     4th order 
      zcoup_e1a=zcoup_e1a + (p(4)*pc(1,n)*zc(3,n) + p(5)*pc(2,n)*zc(2,n)
     & + p(6)*pc(3,n)*zc(1,n))

!     5th order 
      zcoup_e1a=zcoup_e1a + (p(7)*pc(4,n)*zc(1,n) + p(8)*pc(3,n)*zc(2,n)
     & +p(9)*pc(2,n)*zc(3,n)+p(10)*pc(1,n)*zc(4,n)
     & +p(11)*pc(1,n)*zc(5,n))

!     6th order 
      zcoup_e1a = zcoup_e1a + 
     &    p(12)*pc(5,n)*zc(1,n)+p(13)*pc(4,n)*zc(2,n)
     &   +p(14)*pc(3,n)*zc(3,n)
     &    + p(15)*pc(2,n)*zc(4,n)+p(16)*pc(2,n)*zc(5,n)
     &    + pc(1,n)*(p(17)*zc(6,n)+p(18)*zc(7,n))


      end      

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between a 
c and b modes up_ to sixth order
c uses not more than 15 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 4
c 6.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Ev_ab(v_aa, n, par_, j, iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien
      

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=15)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v_aa !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo
!.....

!..   second order
      v_aa = p(1) * pa(1,n) * pb(1,n)

      !third order
      v_aa =v_aa + p(2) * pa(2,n) * pb(1,n) + p(3) * pa(1,n) * pb(2,n)

!..   fourth order
      v_aa =v_aa + p(4) * pa(3,n) * pb(1,n) + p(5) * pa(2,n) * pb(2,n)
     &     + p(6) * pa(1,n) * pb(3,n)

!..   fifth order
      v_aa =v_aa + p(7) * pa(4,n) * pb(1,n) + p(8) * pa(3,n) * pb(2,n)
     &     + p(9) * pa(2,n) * pb(3,n) + p(10) * pa(1,n) * pb(4,n)

!..   sixths order
      v_aa =v_aa + p(11) * pa(5,n) * pb(1,n) + p(12) * pa(4,n) * pb(2,n)
     &     + p(13) * pa(3,n) * pb(3,n) + p(14) * pa(2,n) * pb(4,n)
     &     + p(15) * pa(1,n) * pb(5,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c a and c modes up_ to sixth order
c uses not more than 15 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 4
c 6.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine Ev_ac(v_aa, n, par_, j, iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=15)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v_aa !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!..   second order
      v_aa = p(1) * pa(1,n) * pc(1,n)

c     third order
      v_aa =v_aa + p(2) * pa(2,n) * pc(1,n) + p(3) * pa(1,n) * pc(2,n)

!..   fourth order
      v_aa =v_aa + p(4) * pa(3,n) * pc(1,n) + p(5) * pa(2,n) * pc(2,n)
     &     + p(6) * pa(1,n) * pc(3,n)

!..   fifth order
      v_aa =v_aa + p(7) * pa(4,n) * pc(1,n) + p(8) * pa(3,n) * pc(2,n)
     &     + p(9) * pa(2,n) * pc(3,n) + p(10) * pa(1,n) * pc(4,n)

!..   sixths order
      v_aa =v_aa + p(11) * pa(5,n) * pc(1,n) + p(12) * pa(4,n) * pc(2,n)
     &     + p(13) * pa(3,n) * pc(3,n) + p(14) * pa(2,n) * pc(4,n)
     &     + p(15) * pa(1,n) * pc(5,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between b
c and c modes up_ to sixth order
c uses not more than 15 parameters
c 2.: 1
c 3.: 2
c 4.: 3
c 5.: 4
c 6.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine Ev_bc(v_aa, n, par_, j, iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=15)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v_aa !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'


c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(1,j)
      if (iref.eq.1) then
        ii=2  !2-mod-kop ???
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!..   second order
      v_aa = p(1) * pb(1,n) * pc(1,n)

c     third order
      v_aa =v_aa + p(2) * pb(2,n) * pc(1,n) + p(3) * pb(1,n) * pc(2,n)

!..   fourth order
      v_aa =v_aa + p(4) * pb(3,n) * pc(1,n) + p(5) * pb(2,n) * pc(2,n)
     &     + p(6) * pb(1,n) * pc(3,n)

!..   fifth order
      v_aa =v_aa + p(7) * pb(4,n) * pc(1,n) + p(8) * pb(3,n) * pc(2,n)
     &     + p(9) * pb(2,n) * pc(3,n) + p(10) * pb(1,n) * pc(4,n)

!..   sixths order
      v_aa =v_aa + p(11) * pb(5,n) * pc(1,n) + p(12) * pb(4,n) * pc(2,n)
     &     + p(13) * pb(3,n) * pc(3,n) + p(14) * pb(2,n) * pc(4,n)
     &     + p(15) * pb(1,n) * pc(5,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c all three e-modes
c uses not more than 7 parameters
c 3.: 1
c 4.: 6
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_veee(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=7)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

c     third order
      v = p(1) * veee(1,n)

!..   fourth order
      v = v + p(2)*veee(2,n) + p(3)*veee(3,n) + p(4)*veee(4,n) 
     &     + p(5)*veee(5,n) + p(6)*veee(6,n) + p(7)*veee(7,n)

      end


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c all three e-modes
c uses not more than 15 parameters
c 3.: 3
c 4.: 12
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupeee(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=15)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

c     third order
      wcoup = p(1) * weee(1,n) + p(2)*weee(2,n) + p(3)*weee(3,n) 

!..   fourth order
      wcoup = wcoup + p(4)*weee(4,n) + p(5)*weee(5,n) + p(6)*weee(6,n)
     &     + p(7)*weee(7,n) + p(8)*weee(8,n) + p(9)*weee(9,n)
     &     + p(10)*weee(10,n) + p(11)*weee(11,n) + p(12)*weee(12,n)
     &     + p(13)*weee(13,n) + p(14)*weee(14,n) + p(15)*weee(15,n)  
      end


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c all three e-modes
c uses not more than 15 parameters
c 3.: 3
c 4.: 12
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupeee(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=15)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

c     third order
      zcoup = p(1) * zeee(1,n) + p(2)*zeee(2,n) + p(3)*zeee(3,n) 

!..   fourth order
      zcoup = zcoup + p(4)*zeee(4,n) + p(5)*zeee(5,n) + p(6)*zeee(6,n)
     &     + p(7)*zeee(7,n) + p(8)*zeee(8,n) + p(9)*zeee(9,n)
     &     + p(10)*zeee(10,n) + p(11)*zeee(11,n) + p(12)*zeee(12,n)
     &     + p(13)*zeee(13,n) + p(14)*zeee(14,n) + p(15)*zeee(15,n)  
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c all three a-modes
c uses not more than 4 parameters
c 3.: 1
c 4.: 3
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vaaa(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=4)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

c     third order
      v = p(1)*pa(1,n)*pb(2,n)*pc(1,n)

!..   fourth order
      v = v + p(2)*(pa(1,n)**2)*pb(2,n)*pc(1,n)
      v = v + p(3)*pa(1,n)*(pb(2,n)**2)*pc(1,n)
      v = v + p(4)*pa(1,n)*pb(2,n)*(pc(1,n)**2)
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c a,ea,eb
c uses not more than 4 parameters
c 3.: 1
c 4.: 3
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vaeaeb(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=4)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     3nd order
      v = p(1)*pa(1,n)*veeab(1,n)

!     4th order
      v = v + p(2)*pa(1,n)*veeab(2,n) + p(3)*pa(1,n)*veeab(3,n)
     &     + p(4)*pa(2,n)*veeab(1,n)
      end


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c a,ea,ec
c uses not more than 4 parameters
c 3.: 1
c 4.: 3
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vaeaec(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=4)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     3nd order
      v = p(1)*pa(1,n)*veeac(1,n)

!     4th order
      v = v + p(2)*pa(1,n)*veeac(2,n) + p(3)*pa(1,n)*veeac(3,n)
     &     + p(4)*pa(2,n)*veeac(1,n)
      end


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between  
c a,eb,ec
c uses not more than 4 parameters
c 3.: 1
c 4.: 3
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vaebec(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=4)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     3nd order
      v = p(1)*pa(1,n)*veebc(1,n)

!     4th order
      v = v + p(2)*pa(1,n)*veebc(2,n) + p(3)*pa(1,n)*veebc(3,n)
     &     + p(4)*pa(2,n)*veebc(1,n)
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c b,ea,eb
c uses not more than 4 parameters
c 3.: 1
c 4.: 3
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vbeaeb(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=4)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     3nd order
      v = p(1)*pb(1,n)*veeab(1,n)

!     4th order
      v = v + p(2)*pb(1,n)*veeab(2,n) + p(3)*pb(1,n)*veeab(3,n)
     &     + p(4)*pb(2,n)*veeab(1,n)
      end


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c b,ea,ec
c uses not more than 4 parameters
c 3.: 1
c 4.: 3
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vbeaec(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=4)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     3nd order
      v = p(1)*pb(1,n)*veeac(1,n)

!     4th order
      v = v + p(2)*pb(1,n)*veeac(2,n) + p(3)*pb(1,n)*veeac(3,n)
     &     + p(4)*pb(2,n)*veeac(1,n)
      end


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between  
c b,eb,ec
c uses not more than 4 parameters
c 3.: 1
c 4.: 3
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vbebec(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=4)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     3nd order
      v = p(1)*pb(1,n)*veebc(1,n)

!     4th order
      v = v + p(2)*pb(1,n)*veebc(2,n) + p(3)*pb(1,n)*veebc(3,n)
     &     + p(4)*pb(2,n)*veebc(1,n)
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c c,ea,eb
c uses not more than 4 parameters
c 3.: 1
c 4.: 3
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vceaeb(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=4)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     3nd order
      v = p(1)*pc(1,n)*veeab(1,n)

!     4th order
      v = v + p(2)*pc(1,n)*veeab(2,n) + p(3)*pc(1,n)*veeab(3,n)
     &     + p(4)*pc(2,n)*veeab(1,n)
      end


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c c,ea,ec
c uses not more than 4 parameters
c 3.: 1
c 4.: 3
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vceaec(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=4)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     3nd order
      v = p(1)*pc(1,n)*veeac(1,n)

!     4th order
      v = v + p(2)*pc(1,n)*veeac(2,n) + p(3)*pc(1,n)*veeac(3,n)
     &     + p(4)*pc(2,n)*veeac(1,n)
      end


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between  
c c,eb,ec
c uses not more than 4 parameters
c 3.: 1
c 4.: 3
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vcebec(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=4)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     3nd order
      v = p(1)*pc(1,n)*veebc(1,n)

!     4th order
      v = v + p(2)*pc(1,n)*veebc(2,n) + p(3)*pc(1,n)*veebc(3,n)
     &     + p(4)*pc(2,n)*veebc(1,n)
      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c a,ea,eb
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupaeaeb(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pa(1,n)*weeab(1,n)

!     4th order
      wcoup = wcoup + p(2)*pa(1,n)*weeab(2,n) + p(3)*pa(1,n)*weeab(3,n)
     &     + p(4)*pa(1,n)*weeab(4,n) + p(5)*pa(1,n)*weeab(5,n)
     &     + p(6)*pa(2,n)*weeab(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c a,ea,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupaeaec(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pa(1,n)*weeac(1,n)

!     4th order
      wcoup = wcoup + p(2)*pa(1,n)*weeac(2,n) + p(3)*pa(1,n)*weeac(3,n)
     &     + p(4)*pa(1,n)*weeac(4,n) + p(5)*pa(1,n)*weeac(5,n)
     &     + p(6)*pa(2,n)*weeac(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c a,eb,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupaebec(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pa(1,n)*weebc(1,n)

!     4th order
      wcoup = wcoup + p(2)*pa(1,n)*weebc(2,n) + p(3)*pa(1,n)*weebc(3,n)
     &     + p(4)*pa(1,n)*weebc(4,n) + p(5)*pa(1,n)*weebc(5,n)
     &     + p(6)*pa(2,n)*weebc(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c b,ea,eb
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupbeaeb(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pb(1,n)*weeab(1,n)

!     4th order
      wcoup = wcoup + p(2)*pb(1,n)*weeab(2,n) + p(3)*pb(1,n)*weeab(3,n)
     &     + p(4)*pb(1,n)*weeab(4,n) + p(5)*pb(1,n)*weeab(5,n)
     &     + p(6)*pb(2,n)*weeab(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c b,ea,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupbeaec(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pb(1,n)*weeac(1,n)

!     4th order
      wcoup = wcoup + p(2)*pb(1,n)*weeac(2,n) + p(3)*pb(1,n)*weeac(3,n)
     &     + p(4)*pb(1,n)*weeac(4,n) + p(5)*pb(1,n)*weeac(5,n)
     &     + p(6)*pb(2,n)*weeac(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c b,eb,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupbebec(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pb(1,n)*weebc(1,n)

!     4th order
      wcoup = wcoup + p(2)*pb(1,n)*weebc(2,n) + p(3)*pb(1,n)*weebc(3,n)
     &     + p(4)*pb(1,n)*weebc(4,n) + p(5)*pb(1,n)*weebc(5,n)
     &     + p(6)*pb(2,n)*weebc(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c c,ea,eb
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupceaeb(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pc(1,n)*weeab(1,n)

!     4th order
      wcoup = wcoup + p(2)*pc(1,n)*weeab(2,n) + p(3)*pc(1,n)*weeab(3,n)
     &     + p(4)*pc(1,n)*weeab(4,n) + p(5)*pc(1,n)*weeab(5,n)
     &     + p(6)*pc(2,n)*weeab(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c c,ea,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupceaec(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pc(1,n)*weeac(1,n)

!     4th order
      wcoup = wcoup + p(2)*pc(1,n)*weeac(2,n) + p(3)*pc(1,n)*weeac(3,n)
     &     + p(4)*pc(1,n)*weeac(4,n) + p(5)*pc(1,n)*weeac(5,n)
     &     + p(6)*pc(2,n)*weeac(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c c,eb,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupcebec(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pc(1,n)*weebc(1,n)

!     4th order
      wcoup = wcoup + p(2)*pc(1,n)*weebc(2,n) + p(3)*pc(1,n)*weebc(3,n)
     &     + p(4)*pc(1,n)*weebc(4,n) + p(5)*pc(1,n)*weebc(5,n)
     &     + p(6)*pc(2,n)*weebc(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c a,ea,eb
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupaeaeb(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pa(1,n)*zeeab(1,n)

!     4th order
      zcoup = zcoup + p(2)*pa(1,n)*zeeab(2,n) + p(3)*pa(1,n)*zeeab(3,n)
     &     + p(4)*pa(1,n)*zeeab(4,n) + p(5)*pa(1,n)*zeeab(5,n)
     &     + p(6)*pa(2,n)*zeeab(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c a,ea,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupaeaec(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pa(1,n)*zeeac(1,n)

!     4th order
      zcoup = zcoup + p(2)*pa(1,n)*zeeac(2,n) + p(3)*pa(1,n)*zeeac(3,n)
     &     + p(4)*pa(1,n)*zeeac(4,n) + p(5)*pa(1,n)*zeeac(5,n)
     &     + p(6)*pa(2,n)*zeeac(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c a,eb,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupaebec(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pa(1,n)*zeebc(1,n)

!     4th order
      zcoup = zcoup + p(2)*pa(1,n)*zeebc(2,n) + p(3)*pa(1,n)*zeebc(3,n)
     &     + p(4)*pa(1,n)*zeebc(4,n) + p(5)*pa(1,n)*zeebc(5,n)
     &     + p(6)*pa(2,n)*zeebc(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c b,ea,eb
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupbeaeb(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pb(1,n)*zeeab(1,n)

!     4th order
      zcoup = zcoup + p(2)*pb(1,n)*zeeab(2,n) + p(3)*pb(1,n)*zeeab(3,n)
     &     + p(4)*pb(1,n)*zeeab(4,n) + p(5)*pb(1,n)*zeeab(5,n)
     &     + p(6)*pb(2,n)*zeeab(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c b,ea,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupbeaec(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pb(1,n)*zeeac(1,n)

!     4th order
      zcoup = zcoup + p(2)*pb(1,n)*zeeac(2,n) + p(3)*pb(1,n)*zeeac(3,n)
     &     + p(4)*pb(1,n)*zeeac(4,n) + p(5)*pb(1,n)*zeeac(5,n)
     &     + p(6)*pb(2,n)*zeeac(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c b,eb,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupbebec(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pb(1,n)*zeebc(1,n)

!     4th order
      zcoup = zcoup + p(2)*pb(1,n)*zeebc(2,n) + p(3)*pb(1,n)*zeebc(3,n)
     &     + p(4)*pb(1,n)*zeebc(4,n) + p(5)*pb(1,n)*zeebc(5,n)
     &     + p(6)*pb(2,n)*zeebc(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c c,ea,eb
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupceaeb(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pc(1,n)*zeeab(1,n)

!     4th order
      zcoup = zcoup + p(2)*pc(1,n)*zeeab(2,n) + p(3)*pc(1,n)*zeeab(3,n)
     &     + p(4)*pc(1,n)*zeeab(4,n) + p(5)*pc(1,n)*zeeab(5,n)
     &     + p(6)*pc(2,n)*zeeab(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c c,ea,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupceaec(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pc(1,n)*zeeac(1,n)

!     4th order
      zcoup = zcoup + p(2)*pc(1,n)*zeeac(2,n) + p(3)*pc(1,n)*zeeac(3,n)
     &     + p(4)*pc(1,n)*zeeac(4,n) + p(5)*pc(1,n)*zeeac(5,n)
     &     + p(6)*pc(2,n)*zeeac(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c c,eb,ec
c uses not more than 6 parameters
c 3.: 1
c 4.: 5
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupcebec(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=6)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pc(1,n)*zeebc(1,n)

!     4th order
      zcoup = zcoup + p(2)*pc(1,n)*zeebc(2,n) + p(3)*pc(1,n)*zeebc(3,n)
     &     + p(4)*pc(1,n)*zeebc(4,n) + p(5)*pc(1,n)*zeebc(5,n)
     &     + p(6)*pc(2,n)*zeebc(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c a,b,ea
c uses not more than 1 parameters
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vabea(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=1)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     4th order
      v = p(1)*pa(1,n)*pb(1,n)*va(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c a,c,ea
c uses not more than 1 parameters
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vacea(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=1)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     4th order
      v = p(1)*pa(1,n)*pc(1,n)*va(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c b,c,ea
c uses not more than 1 parameters
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vbcea(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=1)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     4th order
      v = p(1)*pb(1,n)*pc(1,n)*va(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c a,b,eb
c uses not more than 1 parameters
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vabeb(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=1)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     4th order
      v = p(1)*pa(1,n)*pb(1,n)*vb(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c a,c,eb
c uses not more than 1 parameters
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vaceb(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=1)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     4th order
      v = p(1)*pa(1,n)*pc(1,n)*vb(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c b,c,eb
c uses not more than 1 parameters
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vbceb(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=1)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     4th order
      v = p(1)*pb(1,n)*pc(1,n)*vb(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c a,b,ec
c uses not more than 1 parameters
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vabec(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=1)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     4th order
      v = p(1)*pa(1,n)*pb(1,n)*vc(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c a,c,ec
c uses not more than 1 parameters
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vacec(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=1)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     4th order
      v = p(1)*pa(1,n)*pc(1,n)*vc(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate V matrix elements for the coupling between 
c b,c,ec
c uses not more than 1 parameters
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine e_vbcec(v,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=1)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision v !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo

!     4th order
      v = p(1)*pb(1,n)*pc(1,n)*vc(1,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c a,b,ea
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupabea(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pa(1,n)*pb(1,n)*wa(1,n)

!     4th order
      wcoup = wcoup + p(2)*pa(1,n)*pb(1,n)*wa(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c a,c,ea
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupacea(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pa(1,n)*pc(1,n)*wa(1,n)

!     4th order
      wcoup = wcoup + p(2)*pa(1,n)*pc(1,n)*wa(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c b,c,ea
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupbcea(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pb(1,n)*pc(1,n)*wa(1,n)

!     4th order
      wcoup = wcoup + p(2)*pb(1,n)*pc(1,n)*wa(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c a,b,eb
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupabeb(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pa(1,n)*pb(1,n)*wb(1,n)

!     4th order
      wcoup = wcoup + p(2)*pa(1,n)*pb(1,n)*wb(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c a,c,eb
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupaceb(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pa(1,n)*pc(1,n)*wb(1,n)

!     4th order
      wcoup = wcoup + p(2)*pa(1,n)*pc(1,n)*wb(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c b,c,eb
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupbceb(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pb(1,n)*pc(1,n)*wb(1,n)

!     4th order
      wcoup = wcoup + p(2)*pb(1,n)*pc(1,n)*wb(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c a,b,ec
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupabec(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pa(1,n)*pb(1,n)*wc(1,n)

!     4th order
      wcoup = wcoup + p(2)*pa(1,n)*pb(1,n)*wc(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c a,c,ec
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupacec(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pa(1,n)*pc(1,n)*wc(1,n)

!     4th order
      wcoup = wcoup + p(2)*pa(1,n)*pc(1,n)*wc(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate W matrix elements for the coupling between 
c b,c,ec
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine wcoupbcec(wcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision wcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      wcoup = p(1)*pb(1,n)*pc(1,n)*wc(1,n)

!     4th order
      wcoup = wcoup + p(2)*pb(1,n)*pc(1,n)*wc(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c a,b,ea
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupabea(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pa(1,n)*pb(1,n)*za(1,n)

!     4th order
      zcoup = zcoup + p(2)*pa(1,n)*pb(1,n)*za(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c a,c,ea
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupacea(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pa(1,n)*pc(1,n)*za(1,n)

!     4th order
      zcoup = zcoup + p(2)*pa(1,n)*pc(1,n)*za(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c b,c,ea
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupbcea(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pb(1,n)*pc(1,n)*za(1,n)

!     4th order
      zcoup = zcoup + p(2)*pb(1,n)*pc(1,n)*za(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c a,b,eb
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupabeb(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pa(1,n)*pb(1,n)*zb(1,n)

!     4th order
      zcoup = zcoup + p(2)*pa(1,n)*pb(1,n)*zb(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c a,c,eb
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupaceb(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pa(1,n)*pc(1,n)*zb(1,n)

!     4th order
      zcoup = zcoup + p(2)*pa(1,n)*pc(1,n)*zb(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c b,c,eb
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupbceb(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pb(1,n)*pc(1,n)*zb(1,n)

!     4th order
      zcoup = zcoup + p(2)*pb(1,n)*pc(1,n)*zb(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c a,b,ec
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupabec(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pa(1,n)*pb(1,n)*zc(1,n)

!     4th order
      zcoup = zcoup + p(2)*pa(1,n)*pb(1,n)*zc(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c a,c,ec
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupacec(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pa(1,n)*pc(1,n)*zc(1,n)

!     4th order
      zcoup = zcoup + p(2)*pa(1,n)*pc(1,n)*zc(2,n)

      end

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c function to generate Z matrix elements for the coupling between 
c b,c,ec
c uses not more than 2 parameters
c 3.: 1
c 4.: 1
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine zcoupbcec(zcoup,n,par_,j,iref)
      implicit none

      integer i !running indecies

      integer n !point

      integer iref, ii, jj !used for refference Hamiltonien

      integer j               !number of given parameters
      integer lnx             !Maximum number of used parameters
      parameter (lnx=2)
      double precision par_(j) !given parameter vector
      double precision p(lnx) !internal parameter vector

      double precision zcoup !Matrix element

      include 'params.incl'
      include 'vwz_MeX.incl'

c     copy parameters to local ones to account for shorter than max.
c     expansions
      do i = 1, lnx
         p(i) = 0.d0
      enddo

c     this takes care of split reference plus correction scheme
      ii=1                
      jj=min(0,j)
      if (iref.eq.1) then
        ii=1
        jj=j
      endif

      do i = ii, jj
         if(j.gt.0) p(i) = par_(i)
      enddo


!     3rd order      
      zcoup = p(1)*pb(1,n)*pc(1,n)*zc(1,n)

!     4th order
      zcoup = zcoup + p(2)*pb(1,n)*pc(1,n)*zc(2,n)

      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE MAKESOMATRIX(NDIAB,SOMATRIX,PST,PARAM,KOORD)
C% T. Weike August 2015
C%
C% Konstruiert die volldimensionale HSO-Matrix aus Parametern und 
C% Koordinaten in C3v und 9 Dimensionen. Das Ungekoppelte System 
C% wird bis zur 4. Ordnung beschrieben. Die Kopplung ist nur in
C% exe und bis zur 2. Ordnung implementiert.
C% ndiab:    Dimension der Matrix (int)
C% soMatrix: Die MAtrix, welche beschrieben Werden soll gemaess dem 
C%           HSO-Modell fuer C3v (double complex(ndiab,ndiab))
C% pst:      Angabe ueber die Parameterblöcke pst(1,x): Startpunkt
C%           und pst(2,x): Laenge des Bereichs (int)
C% param:    Parameter Vektor, enthaelt alle reelen Parameter (double)
C% koord:    9 Koordinaten(x,y,z) bzw (a1,ex,ey) (double)
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine makeSOMatrix(ndiab,    !integer
     f                        soMatrix, !double complex(ndiab,ndiab)
     f                        pst,      !integer(2,*)
     f                        param,    !double precision(*)
     f                        koord)    !double precision(dim,*)
      implicit none
      
c     Konstanten
      integer    p      !Maximale Laenge des Parametervektors
      parameter (p=400)

c     Variablen
      integer ndiab   !Rang der Matrix

c     Es wird die SO-Matrix zurückgegeben
      double complex, dimension(ndiab,ndiab) :: soMatrix

c     Parameter und Koordinaten
      integer,          dimension(2,22) :: pst   !Struktur des param-Vek.
      double precision, dimension(p)    :: param !Vek. der Parameter (Vorfaktoren)
      double precision, dimension(9,1)  :: koord !Koordinaten des System (Ortsvek.)
      double precision, dimension(6,1)  :: hilf  !Hilfsarray fuer Koordinaten
      
c     Laufindex
      integer i
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Initialisierung der SO-Matrix
      soMatrix=0.d0

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Erstellt die 0. Ordnung
       call makeNullteOrdnung(ndiab,soMatrix,pst,param)

c     Erstellt die Ungekoppelten Ordnungen 1-4
c     Este mode (x1,y1,z1)
       call makeSOMatrixDreiD(ndiab,soMatrix,
     f                       pst(1:2, 2: 5),param,koord(1:3,1))
c     Zweite mode (x2,y2,z2)
      call makeSOMatrixDreiD(ndiab,soMatrix,
     f                       pst(1:2, 6: 9),param,koord(4:6,1))
c     Dritte mode (x3,y3,z3)
      call makeSOMatrixDreiD(ndiab,soMatrix,
     f                       pst(1:2,10:13),param,koord(7:9,1))
!      call printComplexMatrix(ndiab,ndiab,soMatrix)

c     Kopplung 1+2
      call makeSOMatrixZweimodenkoplung(ndiab,soMatrix,pst(1:2,14:16),
     f                                  param,koord(1:6,1))

c     Kopplung 1+3
      hilf(1,1) = koord(1,1)
      hilf(2,1) = koord(2,1)
      hilf(3,1) = koord(3,1)
      hilf(4,1) = koord(7,1)
      hilf(5,1) = koord(8,1)
      hilf(6,1) = koord(9,1)
      call makeSOMatrixZweimodenkoplung(ndiab,soMatrix,pst(1:2,17:19),
     f                                  param,hilf)

c     Kopplung 2+3
      call makeSOMatrixZweimodenkoplung(ndiab,soMatrix,pst(1:2,20:22),
     f                                  param,koord(4:9,1))

      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Konstruiert 1. bis 4. Ordnung der HSO-Matrix aus Parametern und 
c     Koordinaten in C3v und 3 Dimensionen (Ungekoppelt)
c     ndiab: Dimension der Matrix
c     soMatrix: Die MAtrix, welche beschrieben Werden soll gemaess dem 
c               HSO-Modell fuer C3v
c     pst:      Angabe ueber die Parameterblöcke pst(1,x): Startpunkt
c               und pst(2,x): Laenge des Bereichs
c     param:    Parameter Vektor, enthaelt alle reelen Parameter
c     koord:    3 Koordinaten(x,y,z) bzw (a1,ex,ey)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine makeSOMatrixDreiD(ndiab,    !integer
     f                             soMatrix, !double complex(ndiab,ndiab)
     f                             pst,      !integer(2,x:y), mit x+4=y
     f                             param,    !double precision(*)
     f                             koord)    !double precision(3,*)
      implicit none
      
c     Parameter
      integer ndiab   !Rang der Matrix

c     Es wird die SO-Matrix zurückgegeben
      double complex, dimension(ndiab,ndiab) :: soMatrix

c     Einheiten Komplexer Zahlen
      double complex im                  !Imaginäre Einheit
      double complex re                  !Reelle Einheit
      parameter (im=dcmplx(0.0d0,1.0d0))
      parameter (re=dcmplx(1.0d0,0.0d0))

c     Parameter und Koordinaten
      integer,          dimension(2 ,4) :: pst    !Struktur des param-Vek.         
      double precision, dimension(400)  :: param  !Vek. der Parameter (Vorfaktoren)
      double precision, dimension(3,1)  :: koord  !Koordinaten des System (Ortsvek.)

c     Komplexe Koordinaten
      double complex m,p                                !1.Ordnung
      double complex zz,mm,pp,rmm,imm,mp,izp,rzp,zp,zm  !2.Ordnung
      double complex mpz,mpp,mpm,rmpp,impp,zzz,zmm,zpp  !3.Ordnung
      double complex rzmm,izmm,izzp,rzzp,zzp,zzm,mmm
      double complex ppp,rppp
      double complex mpzz, mpzp, mpzm,rmpzp,impzp, mppp !4.Ordnung
      double complex mpmm,rmpmm,impmm, mpmp, zzzz, zzzp
      double complex zzzm,rzzzp,izzzp, zzpp, zzmm,rzzmm
      double complex izzmm, zppp, zmmm,rzppp,pppp, mmmm
      double complex rpppp,ipppp                       

c     Laufindices
      integer i,j

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Koordinaten (vermeidet zusaetzliche Rechnungen)
c     1.Ordnung
      p=koord(1,1)*re+koord(2,1)*im
      m=conjg(p)
c     2.Ordnung
       zz=koord(3,1)**2
       zp=koord(3,1)*p
       zm=koord(3,1)*m
      rzp= dble(zp)
      izp=aimag(zp)
       pp=p**2
       mm=m**2
      rmm= dble(mm)
      imm=aimag(mm)
       mp=m*p
c     3.Ordnung
       mpz=mp*koord(3,1)
       mpp=mp*p
       mpm=mp*m
      rmpp=mp*koord(1,1)
      impp=mp*koord(2,1)
       zzz=koord(3,1)*zz
       zzp=koord(3,1)*zp
       zzm=koord(3,1)*zm
      rzzp= dble(zzp)
      izzp=aimag(zzp)
       zpp=koord(3,1)*pp
       zmm=koord(3,1)*mm
      rzmm= dble(zmm)
      izmm=aimag(zmm)
       ppp=p*pp
       mmm=m*mm
      rppp= dble(ppp)
c     4.Ordnung
       mpzz=mp* zz
       mpzp=mp* zp
       mpzm=mp* zm
      rmpzp=mp*rzp
      impzp=mp*izp
       mppp=mp* pp
       mpmm=mp* mm
      rmpmm=mp*rmm
      impmm=mp*imm
       mpmp=mp* mp
       zzzz=koord(3,1)* zzz
       zzzp=koord(3,1)* zzp
       zzzm=koord(3,1)* zzm
      rzzzp=koord(3,1)*rzzp
      izzzp=koord(3,1)*izzp
       zzpp=koord(3,1)* zpp
       zzmm=koord(3,1)* zmm
      rzzmm=koord(3,1)*rzmm
      izzmm=koord(3,1)*izmm
       zppp=koord(3,1)* ppp
       zmmm=koord(3,1)* mmm
      rzppp=koord(3,1)*rppp
       pppp=p*ppp
       mmmm=m*mmm
      rpppp= dble(pppp)
      ipppp=aimag(pppp)
      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     1. Ordnung
      if(pst(2,1).lt.5) return
c     Epsilon(1,0,0,0)
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,1)  )*im * koord(3,1)
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,1)  )*im * koord(3,1)
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,1)  )*im * koord(3,1)
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,1)  )*im * koord(3,1)
c     Iota(1,0,0,0)
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+1)*re * koord(3,1)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,1)+1)*im * koord(3,1)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,1)+1)*im * koord(3,1)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+1)*re * koord(3,1)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+1)*re * koord(3,1)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,1)+1)*im * koord(3,1)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,1)+1)*im * koord(3,1)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+1)*re * koord(3,1)
c     Zeta(0,0,1,0)
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,1)+2)*im * koord(2,1)
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,1)+2)*im * koord(1,1)
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,1)+2)*im * koord(1,1)
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,1)+2)*im * koord(2,1)
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,1)+2)*im * koord(2,1)
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,1)+2)*im * koord(1,1)
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,1)+2)*im * koord(1,1)
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,1)+2)*im * koord(2,1)
c     Theta(0,0,1,0)
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,1)+3)*im * p
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,1)+3)*im * p
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,1)+3)*im * m
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,1)+3)*im * m
c     Eta(0,0,0,1)
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+4)*re * m
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,1)+4)*im * m
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,1)+4)*im * m
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+4)*re * m
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+4)*re * p
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,1)+4)*im * p
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,1)+4)*im * p
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+4)*re * p

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     2. Ordnung
      if(pst(2,2).lt.10) return
c     Epsilon(0,1,0,0) 
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,2)  )*im *  mp
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,2)  )*im *  mp
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,2)  )*im *  mp
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,2)  )*im *  mp
c     Iota(0,1,0,0)                                       
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+1)*re *  mp
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,2)+1)*im *  mp
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,2)+1)*im *  mp
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+1)*re *  mp
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+1)*re *  mp
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,2)+1)*im *  mp
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,2)+1)*im *  mp
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+1)*re *  mp
c     Epsilon(2,0,0,0)                                    
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,2)+2)*im *  zz
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,2)+2)*im *  zz
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,2)+2)*im *  zz
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,2)+2)*im *  zz
c     Iota(2,0,0,0)                                       
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+3)*re *  zz
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,2)+3)*im *  zz
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,2)+3)*im *  zz
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+3)*re *  zz
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+3)*re *  zz
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,2)+3)*im *  zz
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,2)+3)*im *  zz
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+3)*re *  zz
c     Zeta(1,0,1,0)
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+4)*im * izp
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+4)*im * rzp
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+4)*im * rzp
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+4)*im * izp
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+4)*im * izp
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+4)*im * rzp
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+4)*im * rzp
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+4)*im * izp
c     Theta(1,0,1,0)
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+5)*im *  zp
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+5)*im *  zp
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+5)*im *  zm
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+5)*im *  zm
c     Eta(1,0,0,1)                                        
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+6)*re *  zm
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+6)*im *  zm
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+6)*im *  zm
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+6)*re *  zm
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+6)*re *  zp
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+6)*im *  zp
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+6)*im *  zp
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+6)*re *  zp
c     Zeta(0,0,0,1)
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+7)*im * imm
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+7)*im * rmm
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+7)*im * rmm
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+7)*im * imm
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+7)*im * imm
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+7)*im * rmm
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+7)*im * rmm
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+7)*im * imm
c     Theta(0,0,0,1)
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+8)*im *  mm
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+8)*im *  mm
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+8)*im *  pp
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+8)*im *  pp
c     Eta(0,0,1,0)                                        
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+9)*re *  pp
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+9)*im *  pp
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+9)*im *  pp
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+9)*re *  pp
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+9)*re *  mm
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+9)*im *  mm
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+9)*im *  mm
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+9)*re *  mm

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     3.Ordnung
      if(pst(2,3).lt.16) return
c     Epsilon(1,1,0,0)
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,3)   )*im *  mpz
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,3)   )*im *  mpz
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,3)   )*im *  mpz
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,3)   )*im *  mpz
c     Iota(1,1,0,0)                                        
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,3) +1)*re *  mpz
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,3) +1)*im *  mpz
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,3) +1)*im *  mpz
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,3) +1)*re *  mpz
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,3) +1)*re *  mpz
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,3) +1)*im *  mpz
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,3) +1)*im *  mpz
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,3) +1)*re *  mpz
c     Zeta(0,1,1,0)                              
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,3) +2)*im * impp
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,3) +2)*im * rmpp
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,3) +2)*im * rmpp
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,3) +2)*im * impp
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,3) +2)*im * impp
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,3) +2)*im * rmpp
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,3) +2)*im * rmpp
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,3) +2)*im * impp
c     Theta(0,1,1,0)                             
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,3) +3)*im *  mpp
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,3) +3)*im *  mpp
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,3) +3)*im *  mpm
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,3) +3)*im *  mpm
c     Eta(0,1,0,1)                                         
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,3) +4)*re *  mpm
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,3) +4)*im *  mpm
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,3) +4)*im *  mpm
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,3) +4)*re *  mpm
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,3) +4)*re *  mpp
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,3) +4)*im *  mpp
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,3) +4)*im *  mpp
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,3) +4)*re *  mpp
c     Epsilon(3,0,0,0)                                     
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,3)+ 5)*im *  zzz
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,3)+ 5)*im *  zzz
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,3)+ 5)*im *  zzz
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,3)+ 5)*im *  zzz
c     Iota(3,0,0,0)                                        
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,3)+ 6)*re *  zzz
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,3)+ 6)*im *  zzz
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,3)+ 6)*im *  zzz
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,3)+ 6)*re *  zzz
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,3)+ 6)*re *  zzz
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,3)+ 6)*im *  zzz
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,3)+ 6)*im *  zzz
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,3)+ 6)*re *  zzz
c     Zeta(2,0,1,0)                               
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,3)+ 7)*im * izzp
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,3)+ 7)*im * rzzp
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,3)+ 7)*im * rzzp
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,3)+ 7)*im * izzp
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,3)+ 7)*im * izzp
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,3)+ 7)*im * rzzp
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,3)+ 7)*im * rzzp
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,3)+ 7)*im * izzp
c     Theta(2,0,1,0)                              
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,3)+ 8)*im *  zzp
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,3)+ 8)*im *  zzp
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,3)+ 8)*im *  zzm
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,3)+ 8)*im *  zzm
c     Eta(2,0,0,1)                                         
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,3)+ 9)*re *  zzm
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,3)+ 9)*im *  zzm
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,3)+ 9)*im *  zzm
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,3)+ 9)*re *  zzm
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,3)+ 9)*re *  zzp
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,3)+ 9)*im *  zzp
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,3)+ 9)*im *  zzp
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,3)+ 9)*re *  zzp
c     Zeta(1,0,0,1)
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,3)+10)*im * izmm
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,3)+10)*im * rzmm
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,3)+10)*im * rzmm
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,3)+10)*im * izmm
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,3)+10)*im * izmm
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,3)+10)*im * rzmm
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,3)+10)*im * rzmm
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,3)+10)*im * izmm
c     Theta(1,0,0,1)
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,3)+11)*im *  zmm
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,3)+11)*im *  zmm
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,3)+11)*im *  zpp
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,3)+11)*im *  zpp
c     Eta(1,0,1,0)                                         
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,3)+12)*re *  zpp
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,3)+12)*im *  zpp
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,3)+12)*im *  zpp
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,3)+12)*re *  zpp
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,3)+12)*re *  zmm
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,3)+12)*im *  zmm
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,3)+12)*im *  zmm
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,3)+12)*re *  zmm
c     Epsilon(0,0,1,0)                            
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,3)+13)*im * rppp
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,3)+13)*im * rppp
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,3)+13)*im * rppp
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,3)+13)*im * rppp
c     Iota(0,0,1,0)
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,3)+14)*re *  ppp
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,3)+14)*im *  ppp
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,3)+14)*im *  ppp
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,3)+14)*re *  ppp
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,3)+14)*re *  mmm
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,3)+14)*im *  mmm
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,3)+14)*im *  mmm
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,3)+14)*re *  mmm
c     Iota(0,0,0,1)                                        
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,3)+15)*re *  mmm
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,3)+15)*im *  mmm
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,3)+15)*im *  mmm
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,3)+15)*re *  mmm
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,3)+15)*re *  ppp
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,3)+15)*im *  ppp
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,3)+15)*im *  ppp
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,3)+15)*re *  ppp

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     4.Ordnung
      if(pst(2,4).lt.24) return
c     Epsilon(0,2,0,0) 
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,4)  )*im *  mpmp
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,4)  )*im *  mpmp
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,4)  )*im *  mpmp
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,4)  )*im *  mpmp
c     Iota(0,2,0,0)                                       
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,4)+1)*re *  mpmp
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,4)+1)*im *  mpmp
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,4)+1)*im *  mpmp
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,4)+1)*re *  mpmp
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,4)+1)*re *  mpmp
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,4)+1)*im *  mpmp
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,4)+1)*im *  mpmp
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,4)+1)*re *  mpmp
c     Epsilon(2,1,0,0)                                    
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,4)+2)*im *  mpzz
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,4)+2)*im *  mpzz
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,4)+2)*im *  mpzz
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,4)+2)*im *  mpzz
c     Iota(2,1,0,0)                                       
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,4)+3)*re *  mpzz
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,4)+3)*im *  mpzz
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,4)+3)*im *  mpzz
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,4)+3)*re *  mpzz
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,4)+3)*re *  mpzz
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,4)+3)*im *  mpzz
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,4)+3)*im *  mpzz
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,4)+3)*re *  mpzz
c     Zeta(1,1,1,0)
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,4)+4)*im * impzp
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,4)+4)*im * rmpzp
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,4)+4)*im * rmpzp
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,4)+4)*im * impzp
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,4)+4)*im * impzp
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,4)+4)*im * rmpzp
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,4)+4)*im * rmpzp
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,4)+4)*im * impzp
c     Theta(1,1,1,0)
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,4)+5)*im *  mpzp
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,4)+5)*im *  mpzp
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,4)+5)*im *  mpzm
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,4)+5)*im *  mpzm
c     Eta(1,1,0,1)                                        
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,4)+6)*re *  mpzm
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,4)+6)*im *  mpzm
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,4)+6)*im *  mpzm
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,4)+6)*re *  mpzm
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,4)+6)*re *  mpzp
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,4)+6)*im *  mpzp
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,4)+6)*im *  mpzp
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,4)+6)*re *  mpzp
c     Zeta(0,1,0,1)
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,4)+7)*im * impmm
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,4)+7)*im * rmpmm
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,4)+7)*im * rmpmm
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,4)+7)*im * impmm
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,4)+7)*im * impmm
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,4)+7)*im * rmpmm
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,4)+7)*im * rmpmm
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,4)+7)*im * impmm
c     Theta(0,1,0,1)
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,4)+8)*im *  mpmm
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,4)+8)*im *  mpmm
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,4)+8)*im *  mppp
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,4)+8)*im *  mppp
c     Eta(0,1,1,0)                                        
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,4)+9)*re *  mppp
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,4)+9)*im *  mppp
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,4)+9)*im *  mppp
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,4)+9)*re *  mppp
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,4)+9)*re *  mpmm
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,4)+9)*im *  mpmm
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,4)+9)*im *  mpmm
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,4)+9)*re *  mpmm
c     Epsilon(4,0,0,0)                                     
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,4)+10)*im *  zzzz
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,4)+10)*im *  zzzz
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,4)+10)*im *  zzzz
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,4)+10)*im *  zzzz
c     Iota(4,0,0,0)                                        
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,4)+11)*re *  zzzz
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,4)+11)*im *  zzzz
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,4)+11)*im *  zzzz
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,4)+11)*re *  zzzz
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,4)+11)*re *  zzzz
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,4)+11)*im *  zzzz
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,4)+11)*im *  zzzz
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,4)+11)*re *  zzzz
c     Zeta(3,0,1,0)                               
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,4)+12)*im * izzzp
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,4)+12)*im * rzzzp
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,4)+12)*im * rzzzp
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,4)+12)*im * izzzp
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,4)+12)*im * izzzp
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,4)+12)*im * rzzzp
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,4)+12)*im * rzzzp
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,4)+12)*im * izzzp
c     Theta(3,0,1,0)                              
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,4)+13)*im *  zzzp
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,4)+13)*im *  zzzp
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,4)+13)*im *  zzzm
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,4)+13)*im *  zzzm
c     Eta(3,0,0,1)                                         
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,4)+14)*re *  zzzm
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,4)+14)*im *  zzzm
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,4)+14)*im *  zzzm
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,4)+14)*re *  zzzm
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,4)+14)*re *  zzzp
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,4)+14)*im *  zzzp
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,4)+14)*im *  zzzp
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,4)+14)*re *  zzzp
c     Zeta(2,0,0,1)
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,4)+15)*im * izzmm
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,4)+15)*im * rzzmm
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,4)+15)*im * rzzmm
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,4)+15)*im * izzmm
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,4)+15)*im * izzmm
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,4)+15)*im * rzzmm
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,4)+15)*im * rzzmm
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,4)+15)*im * izzmm
c     Theta(2,0,0,1)
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,4)+16)*im *  zzmm
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,4)+16)*im *  zzmm
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,4)+16)*im *  zzpp
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,4)+16)*im *  zzpp
c     Eta(2,0,1,0)                                         
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,4)+17)*re *  zzpp
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,4)+17)*im *  zzpp
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,4)+17)*im *  zzpp
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,4)+17)*re *  zzpp
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,4)+17)*re *  zzmm
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,4)+17)*im *  zzmm
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,4)+17)*im *  zzmm
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,4)+17)*re *  zzmm
c     Epsilon(1,0,1,0)                            
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,4)+18)*im * rzppp
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,4)+18)*im * rzppp
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,4)+18)*im * rzppp
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,4)+18)*im * rzppp
c     Iota(1,0,1,0)
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,4)+19)*re *  zppp
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,4)+19)*im *  zppp
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,4)+19)*im *  zppp
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,4)+19)*re *  zppp
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,4)+19)*re *  zmmm
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,4)+19)*im *  zmmm
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,4)+19)*im *  zmmm
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,4)+19)*re *  zmmm
c     Iota(1,0,0,1)                                        
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,4)+20)*re *  zmmm
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,4)+20)*im *  zmmm
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,4)+20)*im *  zmmm
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,4)+20)*re *  zmmm
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,4)+20)*re *  zppp
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,4)+20)*im *  zppp
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,4)+20)*im *  zppp
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,4)+20)*re *  zppp
c     Zeta(0,0,2,0)
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,4)+21)*im * ipppp
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,4)+21)*im * rpppp
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,4)+21)*im * rpppp
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,4)+21)*im * ipppp
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,4)+21)*im * ipppp
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,4)+21)*im * rpppp
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,4)+21)*im * rpppp
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,4)+21)*im * ipppp
c     Theta(0,0,1,0)
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,4)+22)*im * pppp
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,4)+22)*im * pppp
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,4)+22)*im * mmmm
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,4)+22)*im * mmmm
c     Eta(0,0,0,1)
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,4)+23)*re * mmmm
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,4)+23)*im * mmmm
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,4)+23)*im * mmmm
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,4)+23)*re * mmmm
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,4)+23)*re * pppp
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,4)+23)*im * pppp
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,4)+23)*im * pppp
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,4)+23)*re * pppp
     
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Erstellt die 0. Ordnung des diabatischen HSO-Modells. Diese wird
c     gesondert erstellt, weil die Koordinatenunabhängig ist.
c     ndiab: Dimension der Matrix
c     soMatrix: Die MAtrix, welche beschrieben Werden soll gemaess dem 
c               HSO-Modell fuer C3v
c     pst:      Angabe ueber die Parameterblöcke pst(1,x): Startpunkt
c               und pst(2,x): Laenge des Bereichs
c     param:    Parameter Vektor, enthaelt alle reelen Parameter
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine makeNullteOrdnung(ndiab,    !integer
     f                             soMatrix, !double complex(ndiab,ndiab)
     f                             pst,      !integer(2,*)
     f                             param)    !double precision(*)
      implicit none
      
c     Parameter
      integer ndiab   !Rang der Matrix

c     Einheiten Komplexer Zahlen
      double complex im                  !Imaginäre Einheit
      double complex re                  !Reelle Einheit
      parameter (im=dcmplx(0.0d0,1.0d0))
      parameter (re=dcmplx(1.0d0,0.0d0))

c     Es wird die SO-Matrix zurückgegeben
      double complex, dimension(ndiab,ndiab) :: soMatrix

c     Parameter und Koordinaten
      integer,          dimension(2,1)  :: pst    !Struktur des param-Vek.         
      double precision, dimension(400)  :: param  !Vek. der Parameter (Vorfaktoren)


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     0. Ordnung
      if(pst(2,1).lt.2) return
c     Epsilon(0,0,0,0) 
      soMatrix(1,2)= param(pst(1,1)  )*im
      soMatrix(2,1)=-param(pst(1,1)  )*im
      soMatrix(5,6)= param(pst(1,1)  )*im
      soMatrix(6,5)=-param(pst(1,1)  )*im
c     Iota(0,0,0,0)
      soMatrix(1,4)= param(pst(1,1)+1)*re
      soMatrix(2,4)= param(pst(1,1)+1)*im
      soMatrix(3,5)=-param(pst(1,1)+1)*im
      soMatrix(3,6)=-param(pst(1,1)+1)*re
      soMatrix(4,1)= param(pst(1,1)+1)*re
      soMatrix(4,2)=-param(pst(1,1)+1)*im
      soMatrix(5,3)= param(pst(1,1)+1)*im
      soMatrix(6,3)=-param(pst(1,1)+1)*re

      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Implementier die Kopplungen zwischen jeweils zwei Moden
c     koord(x1,y1,y1;x2,y2,y3): Koordinaten der Beiden MOden (eigentlich
c                               4 Moden da (a+e) kopplung (a+e)
c     soMatrix: Die (ndiab x ndiab)-Matrix bei der die Kopplungsthereme 
c               hinzugefügt werden sollen
c     param: Vektor mit den Vorfaktoren
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine makeSOMatrixZweimodenkoplung(ndiab,    !integer
     f                                        soMatrix, !double complex(ndiab,ndiab)
     f                                        pst,      !integer(2,x:y), mit x+3=y
     f                                        param,    !double precision(*)
     f                                        koord)    !double precision(6,*)
      implicit none
      
c     Parameter
      integer ndiab   !Rang der Matrix

c     Es wird die SO-Matrix zurückgegeben
      double complex, dimension(ndiab,ndiab) :: soMatrix

c     Einheiten Komplexer Zahlen
      double complex im                  !Imaginäre Einheit
      double complex re                  !Reelle Einheit
      parameter (im=dcmplx(0.0d0,1.0d0))
      parameter (re=dcmplx(1.0d0,0.0d0))

c     Parameter und Koordinaten
      integer,          dimension(2,3) :: pst    !Struktur des param-Vek.         
      double precision, dimension(400) :: param  !Vek. der Parameter (Vorfaktoren)
      double precision, dimension(6,1) :: koord  !Koordinaten des System (Ortsvek.)

c     komplexe Koordinaten
      double complex   pa,pb,ma,mb !q+,q-
      double precision za,zb       !qz
      double complex, dimension(9)  :: oZwei !Zweite Ordnung
      double complex, dimension(33) :: oDrei !Dritte Ordnung

      

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Aufbau der Komplexen Koordinaten
      pa=koord(1,1)*re+koord(2,1)*im
      pb=koord(4,1)*re+koord(5,1)*im
      ma=koord(1,1)*re-koord(2,1)*im
      mb=koord(4,1)*re-koord(5,1)*im
      za=koord(3,1)
      zb=koord(6,1)
c     2.Ordnung
      oZwei(1) = pa*mb
      oZwei(2) = pb*ma
      oZwei(3) = za*zb
      oZwei(4) = za*pb
      oZwei(5) = zb*pa
      oZwei(6) = za*mb
      oZwei(7) = zb*ma
      oZwei(8) = pa*pb
      oZwei(9) = ma*mb
c     3.Ordnung
      oDrei( 1) = za*za*zb
      oDrei( 2) = za*zb*zb
      oDrei( 3) = za*za*pb
      oDrei( 4) = za*zb*pa
      oDrei( 5) = za*zb*pb
      oDrei( 6) = za*za*mb
      oDrei( 7) = za*zb*ma
      oDrei( 8) = za*zb*mb
      oDrei( 9) = za*pa*pb
      oDrei(10) = za*pa*mb
      oDrei(11) = za*ma*pb
      oDrei(12) = za*ma*mb
      oDrei(13) = zb*pa*pa
      oDrei(14) = zb*pa*ma
      oDrei(15) = zb*ma*ma
      oDrei(16) = zb*pa*pb
      oDrei(17) = zb*pa*mb
      oDrei(18) = zb*ma*pb
      oDrei(19) = zb*ma*mb
      oDrei(20) = pa*pa*pb
      oDrei(21) = pa*pa*mb
      oDrei(22) = pa*ma*pb
      oDrei(23) = pa*ma*mb
      oDrei(24) = ma*ma*pb
      oDrei(25) = ma*ma*mb
      oDrei(26) = pa*pb*pb
      oDrei(27) = pa*pb*mb
      oDrei(28) = pa*mb*mb
      oDrei(29) = ma*pb*pb
      oDrei(30) = ma*mb*pb
      oDrei(31) = ma*mb*mb
      oDrei(32) = za*mb*mb
      oDrei(33) = za*pb*pb

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     Aufbau der SO-Matrix
c     2. Ordnung
      if(pst(2,1).lt.15) return
c     Epsilon pa*mb
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,1)   )*im* dble(oZwei(1))
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,1)   )*im* dble(oZwei(1))
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,1)   )*im* dble(oZwei(1))
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,1)   )*im* dble(oZwei(1))
c     Iota pa*mb                                      
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+ 1)*re* oZwei(1)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,1)+ 1)*im* oZwei(1)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,1)+ 1)*im* oZwei(1)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+ 1)*re* oZwei(1)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+ 1)*re* oZwei(2)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,1)+ 1)*im* oZwei(2)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,1)+ 1)*im* oZwei(2)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+ 1)*re* oZwei(2)
c     Iota pb*ma                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+ 3)*re* oZwei(2)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,1)+ 3)*im* oZwei(2)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,1)+ 3)*im* oZwei(2)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+ 3)*re* oZwei(2)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+ 3)*re* oZwei(1)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,1)+ 3)*im* oZwei(1)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,1)+ 3)*im* oZwei(1)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+ 3)*re* oZwei(1)
c     Epsilon zb*za                               
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,1)+ 4)*im* oZwei(3)
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,1)+ 4)*im* oZwei(3)
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,1)+ 4)*im* oZwei(3)
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,1)+ 4)*im* oZwei(3)
c     Iota zb*za                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+ 5)*re* oZwei(3)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,1)+ 5)*im* oZwei(3)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,1)+ 5)*im* oZwei(3)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+ 5)*re* oZwei(3)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+ 5)*re* oZwei(3)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,1)+ 5)*im* oZwei(3)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,1)+ 5)*im* oZwei(3)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+ 5)*re* oZwei(3)
c     Zeta za*pb
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,1)+ 6)*im* aimag(oZwei(4))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,1)+ 6)*im*  dble(oZwei(4))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,1)+ 6)*im*  dble(oZwei(4))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,1)+ 6)*im* aimag(oZwei(4))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,1)+ 6)*im* aimag(oZwei(4))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,1)+ 6)*im*  dble(oZwei(4))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,1)+ 6)*im*  dble(oZwei(4))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,1)+ 6)*im* aimag(oZwei(4))
c     Theta za*pb                                 
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,1)+ 7)*im* oZwei(4)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,1)+ 7)*im* oZwei(4)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,1)+ 7)*im* oZwei(6)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,1)+ 7)*im* oZwei(6)
c     Eta za*mb                                       
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+ 8)*re* oZwei(6)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,1)+ 8)*im* oZwei(6)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,1)+ 8)*im* oZwei(6)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+ 8)*re* oZwei(6)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+ 8)*re* oZwei(4)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,1)+ 8)*im* oZwei(4)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,1)+ 8)*im* oZwei(4)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+ 8)*re* oZwei(4)
c     Zeta zb*pa                                  
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,1)+ 9)*im* aimag(oZwei(5))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,1)+ 9)*im*  dble(oZwei(5))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,1)+ 9)*im*  dble(oZwei(5))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,1)+ 9)*im* aimag(oZwei(5))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,1)+ 9)*im* aimag(oZwei(5))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,1)+ 9)*im*  dble(oZwei(5))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,1)+ 9)*im*  dble(oZwei(5))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,1)+ 9)*im* aimag(oZwei(5))
c     Theta zb*pa                                 
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,1)+10)*im* oZwei(5)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,1)+10)*im* oZwei(5)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,1)+10)*im* oZwei(7)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,1)+10)*im* oZwei(7)
c     Eta zb*ma                                       
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+11)*re* oZwei(7)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,1)+11)*im* oZwei(7)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,1)+11)*im* oZwei(7)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+11)*re* oZwei(7)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+11)*re* oZwei(5)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,1)+11)*im* oZwei(5)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,1)+11)*im* oZwei(5)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+11)*re* oZwei(5)
c     Zeta ma*mb
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,1)+12)*im* aimag(oZwei(9))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,1)+12)*im*  dble(oZwei(9))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,1)+12)*im*  dble(oZwei(9))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,1)+12)*im* aimag(oZwei(9))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,1)+12)*im* aimag(oZwei(9))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,1)+12)*im*  dble(oZwei(9))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,1)+12)*im*  dble(oZwei(9))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,1)+12)*im* aimag(oZwei(9))
c     Theta ma*mb
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,1)+13)*im*  oZwei(9)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,1)+13)*im*  oZwei(9)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,1)+13)*im*  oZwei(8)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,1)+13)*im*  oZwei(8)
c     Eta pa*pb                                        
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+14)*re*  oZwei(8)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,1)+14)*im*  oZwei(8)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,1)+14)*im*  oZwei(8)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+14)*re*  oZwei(8)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+14)*re*  oZwei(9)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,1)+14)*im*  oZwei(9)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,1)+14)*im*  oZwei(9)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+14)*re*  oZwei(9)

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c     3. Ordnung
      if(pst(2,2).lt.48) return
c     Epsilon za*za*zb
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,2)   )*im*oDrei(1)
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,2)   )*im*oDrei(1)
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,2)   )*im*oDrei(1)
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,2)   )*im*oDrei(1)
c     Iota za*za*zb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+ 1)*re*oDrei(1)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,2)+ 1)*im*oDrei(1)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,2)+ 1)*im*oDrei(1)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+ 1)*re*oDrei(1)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+ 1)*re*oDrei(1)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,2)+ 1)*im*oDrei(1)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,2)+ 1)*im*oDrei(1)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+ 1)*re*oDrei(1)
c     Epsilon za*zb*zb
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,2)+ 2)*im*oDrei(2)
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,2)+ 2)*im*oDrei(2)
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,2)+ 2)*im*oDrei(2)
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,2)+ 2)*im*oDrei(2)
c     Iota za*zb*zb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+ 3)*re*oDrei(2)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,2)+ 3)*im*oDrei(2)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,2)+ 3)*im*oDrei(2)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+ 3)*re*oDrei(2)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+ 3)*re*oDrei(2)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,2)+ 3)*im*oDrei(2)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,2)+ 3)*im*oDrei(2)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+ 3)*re*oDrei(2) 
c     Zeta za*za*pb
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+ 4)*im*aimag(oDrei(3))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+ 4)*im* dble(oDrei(3))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+ 4)*im* dble(oDrei(3))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+ 4)*im*aimag(oDrei(3))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+ 4)*im*aimag(oDrei(3))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+ 4)*im* dble(oDrei(3))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+ 4)*im* dble(oDrei(3))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+ 4)*im*aimag(oDrei(3))
c     Theta za*za*pb                                 
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+ 5)*im*oDrei(3)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+ 5)*im*oDrei(3)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+ 5)*im*oDrei(6)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+ 5)*im*oDrei(6)
c     Eta za*za*mb                                      
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+ 6)*re*oDrei(6)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+ 6)*im*oDrei(6)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+ 6)*im*oDrei(6)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+ 6)*re*oDrei(6)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+ 6)*re*oDrei(3)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+ 6)*im*oDrei(3)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+ 6)*im*oDrei(3)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+ 6)*re*oDrei(3) 
c     Zeta za*zb*pa
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+ 7)*im*aimag(oDrei(4))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+ 7)*im* dble(oDrei(4))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+ 7)*im* dble(oDrei(4))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+ 7)*im*aimag(oDrei(4))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+ 7)*im*aimag(oDrei(4))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+ 7)*im* dble(oDrei(4))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+ 7)*im* dble(oDrei(4))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+ 7)*im*aimag(oDrei(4))
c     Theta za*zb*pa                                 
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+ 8)*im*oDrei(4)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+ 8)*im*oDrei(4)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+ 8)*im*oDrei(7)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+ 8)*im*oDrei(7)
c     Eta za*zb*ma                                      
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+ 9)*re*oDrei(7)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+ 9)*im*oDrei(7)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+ 9)*im*oDrei(7)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+ 9)*re*oDrei(7)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+ 9)*re*oDrei(4)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+ 9)*im*oDrei(4)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+ 9)*im*oDrei(4)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+ 9)*re*oDrei(4)
c     Zeta za*zb*pb
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+10)*im*aimag(oDrei(5))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+10)*im* dble(oDrei(5))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+10)*im* dble(oDrei(5))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+10)*im*aimag(oDrei(5))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+10)*im*aimag(oDrei(5))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+10)*im* dble(oDrei(5))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+10)*im* dble(oDrei(5))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+10)*im*aimag(oDrei(5))
c     Theta za*zb*pb                                 
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+11)*im*oDrei(5)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+11)*im*oDrei(5)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+11)*im*oDrei(8)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+11)*im*oDrei(8)
c     Eta za*zb*mb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+12)*re*oDrei(8)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+12)*im*oDrei(8)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+12)*im*oDrei(8)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+12)*re*oDrei(8)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+12)*re*oDrei(5)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+12)*im*oDrei(5)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+12)*im*oDrei(5)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+12)*re*oDrei(5)
c     Zeta za*ma*mb
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+13)*im*aimag(oDrei(12))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+13)*im* dble(oDrei(12))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+13)*im* dble(oDrei(12))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+13)*im*aimag(oDrei(12))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+13)*im*aimag(oDrei(12))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+13)*im* dble(oDrei(12))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+13)*im* dble(oDrei(12))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+13)*im*aimag(oDrei(12))
c     Theta za*ma*mb                                 
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+14)*im*oDrei(12)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+14)*im*oDrei(12)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+14)*im*oDrei( 9)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+14)*im*oDrei( 9)
c     Eta za*pa*pb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+15)*re*oDrei( 9)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+15)*im*oDrei( 9)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+15)*im*oDrei( 9)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+15)*re*oDrei( 9)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+15)*re*oDrei(12)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+15)*im*oDrei(12)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+15)*im*oDrei(12)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+15)*re*oDrei(12)
c     Epsilon za*pa*mb
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,2)+16)*im*dble(oDrei(10))
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,2)+16)*im*dble(oDrei(10))
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,2)+16)*im*dble(oDrei(10))
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,2)+16)*im*dble(oDrei(10))
c     Iota za*pa*mb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+17)*re*oDrei(10)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,2)+17)*im*oDrei(10)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,2)+17)*im*oDrei(10)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+17)*re*oDrei(10)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+17)*re*oDrei(11)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,2)+17)*im*oDrei(11)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,2)+17)*im*oDrei(11)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+17)*re*oDrei(11)
c     Iota Za*pb*ma                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+18)*re*oDrei(11)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,1)+18)*im*oDrei(11)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,1)+18)*im*oDrei(11)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+18)*re*oDrei(11)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+18)*re*oDrei(10)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,1)+18)*im*oDrei(10)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,1)+18)*im*oDrei(10)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+18)*re*oDrei(10)
c     Zeta zb*ma*ma
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+19)*im*aimag(oDrei(15))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+19)*im* dble(oDrei(15))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+19)*im* dble(oDrei(15))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+19)*im*aimag(oDrei(15))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+19)*im*aimag(oDrei(15))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+19)*im* dble(oDrei(15))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+19)*im* dble(oDrei(15))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+19)*im*aimag(oDrei(15))
c     Theta zb*ma*ma                                 
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+20)*im*oDrei(15)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+20)*im*oDrei(15)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+20)*im*oDrei(13)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+20)*im*oDrei(13)
c     Eta zb*pa*pa                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+21)*re*oDrei(13)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+21)*im*oDrei(13)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+21)*im*oDrei(13)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+21)*re*oDrei(13)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+21)*re*oDrei(15)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+21)*im*oDrei(15)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+21)*im*oDrei(15)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+21)*re*oDrei(15)
c     Epsilon zb*pa*ma
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,2)+22)*im*oDrei(14)
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,2)+22)*im*oDrei(14)
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,2)+22)*im*oDrei(14)
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,2)+22)*im*oDrei(14)
c     Iota zb*pa*ma                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+23)*re*oDrei(14)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,2)+23)*im*oDrei(14)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,2)+23)*im*oDrei(14)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+23)*re*oDrei(14)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+23)*re*oDrei(14)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,2)+23)*im*oDrei(14)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,2)+23)*im*oDrei(14)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+23)*re*oDrei(14)
c     Zeta zb*ma*mb
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+24)*im*aimag(oDrei(19))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+24)*im* dble(oDrei(19))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+24)*im* dble(oDrei(19))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+24)*im*aimag(oDrei(19))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+24)*im*aimag(oDrei(19))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+24)*im* dble(oDrei(19))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+24)*im* dble(oDrei(19))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+24)*im*aimag(oDrei(19))
c     Theta zb*ma*mb                                 
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+25)*im*oDrei(19)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+25)*im*oDrei(19)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+25)*im*oDrei(16)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+25)*im*oDrei(16)
c     Eta zb*pa*pb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+26)*re*oDrei(16)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+26)*im*oDrei(16)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+26)*im*oDrei(16)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+26)*re*oDrei(16)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+26)*re*oDrei(19)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+26)*im*oDrei(19)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+26)*im*oDrei(19)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+26)*re*oDrei(19)
c     Epsilon zb*pa*mb
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,2)+27)*im*dble(oDrei(17))
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,2)+27)*im*dble(oDrei(17))
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,2)+27)*im*dble(oDrei(17))
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,2)+27)*im*dble(oDrei(17))
c     Iota zb*pa*mb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+28)*re*oDrei(17)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,2)+28)*im*oDrei(17)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,2)+28)*im*oDrei(17)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+28)*re*oDrei(17)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+28)*re*oDrei(18)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,2)+28)*im*oDrei(18)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,2)+28)*im*oDrei(18)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+28)*re*oDrei(18)
c     Iota zb*pb*ma                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+29)*re*oDrei(18)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,1)+29)*im*oDrei(18)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,1)+29)*im*oDrei(18)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+29)*re*oDrei(18)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+29)*re*oDrei(17)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,1)+29)*im*oDrei(17)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,1)+29)*im*oDrei(17)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+29)*re*oDrei(17)
c     Epsilon pa*pa*pb
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,2)+30)*im*dble(oDrei(20))
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,2)+30)*im*dble(oDrei(20))
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,2)+30)*im*dble(oDrei(20))
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,2)+30)*im*dble(oDrei(20))
c     Iota pa*pa*pb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+31)*re*oDrei(20)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,2)+31)*im*oDrei(20)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,2)+31)*im*oDrei(20)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+31)*re*oDrei(20)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+31)*re*oDrei(25)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,2)+31)*im*oDrei(25)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,2)+31)*im*oDrei(25)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+31)*re*oDrei(25)
c     Iota ma*ma*mb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+32)*re*oDrei(25)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,1)+32)*im*oDrei(25)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,1)+32)*im*oDrei(25)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+32)*re*oDrei(25)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+32)*re*oDrei(20)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,1)+32)*im*oDrei(20)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,1)+32)*im*oDrei(20)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+32)*re*oDrei(20)
c     Zeta pa*pa*mb
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+33)*im*aimag(oDrei(21))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+33)*im* dble(oDrei(21))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+33)*im* dble(oDrei(21))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+33)*im*aimag(oDrei(21))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+33)*im*aimag(oDrei(21))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+33)*im* dble(oDrei(21))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+33)*im* dble(oDrei(21))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+33)*im*aimag(oDrei(21))
c     Theta pa*pa*mb                              
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+34)*im*oDrei(21)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+34)*im*oDrei(21)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+34)*im*oDrei(24)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+34)*im*oDrei(24)
c     Eta ma*ma*pb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+35)*re*oDrei(24)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+35)*im*oDrei(24)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+35)*im*oDrei(24)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+35)*re*oDrei(24)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+35)*re*oDrei(21)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+35)*im*oDrei(21)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+35)*im*oDrei(21)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+35)*re*oDrei(21)
c     Zeta pa*ma*pb
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+36)*im*aimag(oDrei(22))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+36)*im* dble(oDrei(22))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+36)*im* dble(oDrei(22))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+36)*im*aimag(oDrei(22))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+36)*im*aimag(oDrei(22))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+36)*im* dble(oDrei(22))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+36)*im* dble(oDrei(22))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+36)*im*aimag(oDrei(22))
c     Theta pa*ma*pb                              
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+37)*im*oDrei(22)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+37)*im*oDrei(22)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+37)*im*oDrei(23)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+37)*im*oDrei(23)
c     Eta ma*pa*mb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+38)*re*oDrei(23)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+38)*im*oDrei(23)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+38)*im*oDrei(23)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+38)*re*oDrei(23)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+38)*re*oDrei(22)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+38)*im*oDrei(22)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+38)*im*oDrei(22)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+38)*re*oDrei(22)
c     Epsilon pa*pb*pb
      soMatrix(1,2)=soMatrix(1,2)+param(pst(1,2)+39)*im*dble(oDrei(26))
      soMatrix(2,1)=soMatrix(2,1)-param(pst(1,2)+39)*im*dble(oDrei(26))
      soMatrix(5,6)=soMatrix(5,6)+param(pst(1,2)+39)*im*dble(oDrei(26))
      soMatrix(6,5)=soMatrix(6,5)-param(pst(1,2)+39)*im*dble(oDrei(26))
c     Iota pa*pb*pb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+40)*re*oDrei(26)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,2)+40)*im*oDrei(26)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,2)+40)*im*oDrei(26)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+40)*re*oDrei(26)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+40)*re*oDrei(31)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,2)+40)*im*oDrei(31)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,2)+40)*im*oDrei(31)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+40)*re*oDrei(31)
c     Iota ma*mb*mb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,1)+41)*re*oDrei(31)
      soMatrix(2,4)=soMatrix(2,4)+param(pst(1,1)+41)*im*oDrei(31)
      soMatrix(3,5)=soMatrix(3,5)-param(pst(1,1)+41)*im*oDrei(31)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,1)+41)*re*oDrei(31)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,1)+41)*re*oDrei(26)
      soMatrix(4,2)=soMatrix(4,2)-param(pst(1,1)+41)*im*oDrei(26)
      soMatrix(5,3)=soMatrix(5,3)+param(pst(1,1)+41)*im*oDrei(26)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,1)+41)*re*oDrei(26)
c     Zeta pa*pb*mb
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+42)*im*aimag(oDrei(27))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+42)*im* dble(oDrei(27))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+42)*im* dble(oDrei(27))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+42)*im*aimag(oDrei(27))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+42)*im*aimag(oDrei(27))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+42)*im* dble(oDrei(27))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+42)*im* dble(oDrei(27))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+42)*im*aimag(oDrei(27))
c     Theta pa*pb*mb                              
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+43)*im*oDrei(27)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+43)*im*oDrei(27)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+43)*im*oDrei(30)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+43)*im*oDrei(30)
c     Eta ma*pb*mb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+44)*re*oDrei(30)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+44)*im*oDrei(30)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+44)*im*oDrei(30)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+44)*re*oDrei(30)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+44)*re*oDrei(27)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+44)*im*oDrei(27)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+44)*im*oDrei(27)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+44)*re*oDrei(27)
c     Zeta pa*mb*mb
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+45)*im*aimag(oDrei(29))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+45)*im* dble(oDrei(29))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+45)*im* dble(oDrei(29))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+45)*im*aimag(oDrei(29))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+45)*im*aimag(oDrei(29))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+45)*im* dble(oDrei(29))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+45)*im* dble(oDrei(29))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+45)*im*aimag(oDrei(29))
c     Theta pa*mb*mb                              
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+46)*im*oDrei(29)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+46)*im*oDrei(29)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+46)*im*oDrei(28)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+46)*im*oDrei(28)
c     Eta ma*pb*pb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+47)*re*oDrei(28)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+47)*im*oDrei(28)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+47)*im*oDrei(28)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+47)*re*oDrei(28)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+47)*re*oDrei(29)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+47)*im*oDrei(29)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+47)*im*oDrei(29)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+47)*re*oDrei(29)
c     Zeta za*mb*mb
      soMatrix(1,3)=soMatrix(1,3)+param(pst(1,2)+45)*im*aimag(oDrei(32))
      soMatrix(2,3)=soMatrix(2,3)-param(pst(1,2)+45)*im* dble(oDrei(32))
      soMatrix(4,5)=soMatrix(4,5)-param(pst(1,2)+45)*im* dble(oDrei(32))
      soMatrix(4,6)=soMatrix(4,6)+param(pst(1,2)+45)*im*aimag(oDrei(32))
      soMatrix(3,1)=soMatrix(3,1)-param(pst(1,2)+45)*im*aimag(oDrei(32))
      soMatrix(3,2)=soMatrix(3,2)+param(pst(1,2)+45)*im* dble(oDrei(32))
      soMatrix(5,4)=soMatrix(5,4)+param(pst(1,2)+45)*im* dble(oDrei(32))
      soMatrix(6,4)=soMatrix(6,4)-param(pst(1,2)+45)*im*aimag(oDrei(32))
c     Theta za*mb*mb                              
      soMatrix(1,5)=soMatrix(1,5)+param(pst(1,2)+46)*im*oDrei(32)
      soMatrix(2,6)=soMatrix(2,6)-param(pst(1,2)+46)*im*oDrei(32)
      soMatrix(5,1)=soMatrix(5,1)-param(pst(1,2)+46)*im*oDrei(33)
      soMatrix(6,2)=soMatrix(6,2)+param(pst(1,2)+46)*im*oDrei(33)
c     Eta za*pb*pb                                     
      soMatrix(1,4)=soMatrix(1,4)+param(pst(1,2)+47)*re*oDrei(33)
      soMatrix(2,4)=soMatrix(2,4)-param(pst(1,2)+47)*im*oDrei(33)
      soMatrix(3,5)=soMatrix(3,5)+param(pst(1,2)+47)*im*oDrei(33)
      soMatrix(3,6)=soMatrix(3,6)-param(pst(1,2)+47)*re*oDrei(33)
      soMatrix(4,1)=soMatrix(4,1)+param(pst(1,2)+47)*re*oDrei(32)
      soMatrix(4,2)=soMatrix(4,2)+param(pst(1,2)+47)*im*oDrei(32)
      soMatrix(5,3)=soMatrix(5,3)-param(pst(1,2)+47)*im*oDrei(32)
      soMatrix(6,3)=soMatrix(6,3)-param(pst(1,2)+47)*re*oDrei(32)

      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE ADIA(VV,N,Q,P,NPAR,PST,SKIP,Y)
C%
C% determines the adiabatic energies by diagonalizing diabatic matrix.
C% The Eingenvalues are sorted according to the best fitting ordering
C% of the CI vectors.
C%
C% ATTENTION: The interface has changed. To sort by the ci's,
C%            the datavalue of the current points are given
C%
C% input variables:
C% n:    number of point (int)
C% q:    coordinates of the point (double[nq])
C% npar: number of parameters (int)
C% pst:  pointers for parameters (int[np])
C% skip: .false. if everything should be done
C%
C% output variables:
C% vv: firtst nstat energies and than nn cis
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine adia(vv,n,q,p,npar,pst,skip,y)
      implicit none

      include 'states.incl'  !some system parameters
      include 'damp.incl' !somthing for fixing black holes

      integer i,j !running indices
      
      integer iref !getting correction or refference

      double precision eref(nstat,nstat)  !refference (first order)     
      double precision ecorr(nstat,nstat) !correction 2. to wanted order
      double precision e(nstat,nstat)     !full diabatic matrix

      double precision e1(nstat,nstat) !work array 
      double precision e2(nstat,nstat) !work array for diab. matrix
      double precision v1(nstat)       !work array for eigenvalues

      integer n                !current point
      double precision q(qn)   !coordinates
      integer npar             !number of parameters
      integer pst(2,np)        !pointers for parameters
      double precision p(npar) !parameters

      double precision u(nstat,nstat)    !ci-vectors
      double precision uref(nstat,nstat) !ci-vectors
      double precision vv(ntot)          !fitted data
      double precision v(nstat)          !eigen values
      double precision dv(nstat)          !eigen values
      double precision vref(nstat)       !eigen values

      double precision y(ntot)  !read data values

      double precision damp !for fixing black holes

c     Parameters for damping function
      double precision minimum
      double precision width
      parameter(minimum = -0.2, width =0.1)

      logical skip

      skip=.false.

      write(6,*)
      write(6,*) 'q'
      do i=1,9
         write(6,*) q(i)
      enddo

c     compute reference model:
      iref=0
      call diab(eref,n,q,p,pst,iref,npar)
!      if (skip) return
      write(6,*)
      write(6,*) 'eref'
      write(6,*) eref

c     compute correction  model:
      iref=1
      call diab(ecorr,n,q,p,pst,iref,npar)
      write(6,*)
      write(6,*) 'ecorr'
      write(6,*) ecorr

c     copy onto working matrices e1 and e2:
      if (skip) return
      do i=1,nstat
        do j=1,nstat
          e(i,j)=eref(i,j)+ecorr(i,j)
        enddo
      enddo

      call rs(nstat,nstat,eref,vref,1,uref)
      call rs(nstat,nstat,e,v,1,u)

      do i=1,3
         dv(i) = v(i) - vref(i)
      enddo

      write(6,*)
      write(6,*) 'difference'
      write(6,'(3f12.5)') dv

      write(6,*)
      write(6,*) 'total'
      write(6,'(3f12.5)') v

      write(6,*)
      write(6,*) 'refference model'
      write(6,'(3f12.5)') vref

c     now compute full model with damped correction matrix:
      damp = ((tanh((dv(1)-minimum)/width) + 1.)/2.)

c      damp=1.
c      if(dv(1)<minimum) damp = 0.
      do i=1,nstat
        do j=1,nstat
          e(i,j)=eref(i,j)+damp*ecorr(i,j)
        enddo
      enddo

      write(6,*)
      write(6,*) 'e'
      write(6,*) e

c     numeric calculation
      v(1)=0.d0
      call rs(nstat,nstat,e,v,1,u)

c      do i=1,3
c      if(v(i) - vref(i) .lt. -1.) then
c      goto 79
c      endif
c      enddo
c      goto 78
c      
c 79   continue
c      do i=1,3
c      v(i) = vref(i)
c      enddo
c      
c 78   continue

      
c      do i=1,3
c         dv(i) = v(i) - vref(i)
c      enddo
c      do i=1,3
c         v(i) = vref(i) + ((tanh((dv(i)-minimum)/width) + 1.)/2.)*dv(i)
c      enddo
      

c     Writing the energies
      do i=1,nstat
         vv(i)=v(i)
      enddo      
 
c     Writing the cis
      do i=1,nstat
         do j=1,nstat
            vv(nstat*i+j)
!     $           = u(j,i)**2
     $           = u(j,i)
         enddo
      enddo

      call correctPhase(vv(nstat+1),nstat)
      call correctPhase(vv(nstat+4),nstat)
      call correctPhase(vv(nstat+7),nstat)
c      call sortCI(vv(nstat+1),y(nstat+1))

      return

      end

!----------------------------------------------------------------------
!  determine appropriate finite differences for each parameter:
      subroutine pdiff(p,dp,npar)
      implicit none
      integer i, npar
      double precision p(npar), dp(npar), d, thr

      d=1.d-3
      thr=1.d-12

      do i=1,npar
        dp(i)=abs(p(i)*d)
	if (dp(i).lt.thr) dp(i)=thr
      enddo

      end

!----------------------------------------------------------------------
!  compute gradient of adiabatic energies with respect to parameters:
      subroutine vgrad(v,dv,n,q,p,dp,npar,ma,pst,skip,y)
      implicit none

      include 'states.incl'

      integer n, i, j, npar
      integer ma(npar), pst(2,np)
      double precision v(ntot), dv(ntot,*), q(qn),p(npar)
      double precision dp(npar)
      double precision y(ntot)
      logical skip

!     !$omp parallel do default(shared) private(i,j,v,p) 

      do i=1,npar
         do j=1,ntot
            dv(j,i)=0.d0
         enddo
         if (abs(ma(i)).eq.1.d0) then
            p(i)=p(i)+dp(i)
            
            call adia(v,n,q,p,npar,pst,skip,y)
            
            if (skip) return ! be careful, parameter might be changed
            do j=1,ntot
               dv(j,i)=v(j)
            enddo
            p(i)=p(i)-2.d0*dp(i)
            call adia(v,n,q,p,npar,pst,skip,y)
            if (skip) return ! be careful, parameter might be changed
            do j=1,ntot
               dv(j,i)=(dv(j,i)-v(j))/(2.d0*dp(i))
            enddo
            p(i)=p(i)+dp(i)
         endif
      enddo

!!$omp end parallel do

      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C% SUBROUTINE ADIASO(VV,N,Q,P,NPAR,PST,SKIP)
C%
C% determines the adiabatic energies by diagonalizing diabatic matrix.
C% The Eingenvalues are sorted according to the best fitting ordering
C% of the CI vectors.
C%
C% ATTENTION: No hybrid fit implemented
C%
C% input variables:
C% n:    number of point (int)
C% q:    coordinates of the point (double[nq])
C% npar: number of parameters (int)
C% pst:  pointers for parameters (int[np])
C% skip: .false. if everything should be done
C%
C% output variables:
C% vv: firtst nstat energies
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine adiaSO(vv,n,q,p,npar,pst,skip)
      implicit none

      include 'states.incl'     !some system parameters
      include 'damp.incl'    !somthing for fixing black holes
      
      integer i,j !running indices
      
      integer iref !getting correction or refference

      double precision eref(nstat,nstat)  !refference (first order)     
      double precision ecorr(nstat,nstat) !correction 2. to wanted order
      double precision e(nstat,nstat)     !full diabatic matrix
      double complex eso(2*nstat,2*nstat) !diabatic so matrix

      double precision e1(nstat,nstat) !work array 
      double precision e2(nstat,nstat) !work array for diab. matrix
      double precision v1(nstat)       !work array for eigenvalues

      integer n                !current point
      double precision q(qn)   !coordinates
      integer npar             !number of parameters
      integer pst(2,np)        !pointers for parameters
      double precision p(npar) !parameters

      double precision vv(2*nstat)      !fitted data
      double precision v(2*nstat)       !eigen values
      double complex u(2*nstat,2*nstat) !transformation matrix

      double precision damp !for fixing black holes

      logical skip

      include 'params.incl'
      include 'vwz_MeX.incl' ! koordinates

      skip=.false.

c     compute reference model:
      iref=0
      call diab(eref,n,q,p,pst,iref,npar)
!      if (skip) return

c     compute correction  model:
      iref=1
      call diab(ecorr,n,q,p,pst,iref,npar)

c     copy onto working matrices e1 and e2:
      if (skip) return
      do i=1,nstat
        do j=1,nstat
          e1(i,j)=eref(i,j)
          e2(i,j)=eref(i,j)+ecorr(i,j)
        enddo
      enddo

      v1(1)=e1(1,1)
      damp=1.d0

c     now compute full model with damped correction matrix:
      do i=1,nstat
        do j=1,nstat
          e(i,j)=eref(i,j)+damp*ecorr(i,j)
        enddo
      enddo

      do i=1,2*nstat
         do j=1,2*nstat
            eso(i,j)=0.d0
         enddo
      enddo

      call makeSOMatrix(2*nstat,eso,pst(1,283),p,
     &     sokoord(1,n))

c     adding first block of vib coup mat.
      eso(1,1) = eso(1,1) + e(1,1)
      eso(2,1) = eso(2,1) + e(2,1)
      eso(3,1) = eso(3,1) + e(3,1)
      eso(1,2) = eso(1,2) + e(1,2)
      eso(2,2) = eso(2,2) + e(2,2)
      eso(3,2) = eso(3,2) + e(3,2)
      eso(1,3) = eso(1,3) + e(1,3)
      eso(2,3) = eso(2,3) + e(2,3)
      eso(3,3) = eso(3,3) + e(3,3)

c     adding second block of vib coup mat.
      eso(4,4) = eso(4,4) + e(3,3)
      eso(5,4) = eso(5,4) + e(3,2)
      eso(6,4) = eso(6,4) + e(3,1)
      eso(4,5) = eso(4,5) + e(2,3)
      eso(5,5) = eso(5,5) + e(2,2)
      eso(6,5) = eso(6,5) + e(2,1)
      eso(4,6) = eso(4,6) + e(1,3)
      eso(5,6) = eso(5,6) + e(1,2)
      eso(6,6) = eso(6,6) + e(1,1)

c     numeric calculation
      v(1)=0.d0
      call hermdiag(6,6,eso,v,0,u)

c     Writing the energies
      do i=1,2*nstat
         vv(i)=v(i)
      enddo

      return

      end
      
!--------------------------------------------------------------------
! evaluate switching function:
      subroutine ch3clswitch(damp,e,e0ref,e0rho)
      implicit none
      real*8 damp  !damping factor in interval [0,1] on output
      real*8 e     !energy on input
      real*8 e0ref !energy where damping sets in
      real*8 e0rho !width over which damping factor goes to zero
      real*8 eint

!..   damping according to absolute energy:
!      damp=e0rho*(e-e0ref)
!      damp=1.d0-(tanh(damp)+1.d0)*0.5d0

      eint=e-e0ref        !shift energy to onset
      eint=eint/e0rho     !scale energy for required switching range

!      write (6,'(3f12.6)') e, eint

      damp=1.d0
      if (eint.lt.0.d0) return

      damp=0.d0
      if (eint.gt.1.d0) return

!      damp=1.d0 - eint**3 * (2.5d0 - 1.5d0*eint**2)
      damp=1.d0 - eint**2 * (2.d0 - eint**2)
!      damp=(1.d0+cos(eint*3.1415926535897931d0))/2.d0

      end

      subroutine sortStates(v,n)
      implicit none

      include 'states.incl'

      integer n             !length of vector
      double precision v(n) !vector, wich should be sorted
      double precision u(n) !intermidiat save

      double precision tmp !saves the biggest coef
      integer bi           !save index of biggest coef

      integer i,j

      write(6,'(12f7.2)') v
      do j=1,nstat-1
         tmp=0.d0
         do i=nstat*j+j,ntot-nstat,nstat
            if(tmp.lt.v(i)) then
               tmp=v(i)
               bi=(i-j)/nstat
            endif
            write(6,'(i2'': ''f7.2'', ''i2)') bi,v(i),i
         enddo
         do i=1,nstat
            u(j) =v(j)
            v(j) =v(bi)
            v(bi)=u(j)
            u(j *nstat+i)=v(j *nstat+i)
            v(j *nstat+i)=v(bi*nstat+i)
            v(bi*nstat+i)=u(j *nstat+i)
         enddo
      enddo

      end subroutine


      subroutine sortCI(v,y)
      implicit none
      
      include 'states.incl'

      double precision v(nn)
      double precision y(nn)

      integer i

      double precision tmp

      if(y(1).gt.y(2)) then
         if(v(1).gt.v(2)) then
            return
         else
            tmp=v(1)
            v(1)=v(2)
            v(2)=tmp

            tmp=v(4)
            v(4)=v(5)
            v(5)=tmp

            tmp=v(7)
            v(7)=v(8)
            v(8)=tmp
         endif
      else
         if(v(2).gt.v(1)) then
            return
         else
            tmp=v(1)
            v(1)=v(2)
            v(2)=tmp

            tmp=v(4)
            v(4)=v(5)
            v(5)=tmp

            tmp=v(7)
            v(7)=v(8)
            v(8)=tmp
         endif
      endif
      
      end subroutine

      subroutine calculatePotential(e,q,p,npar,pst,skip)
      implicit none

      include 'states.incl'  !some system parameters
      include 'damp.incl' !somthing for fixing black holes

      integer i,j !running indices
      
      integer iref !getting correction or refference

      integer n

      double precision eref(nstat,nstat)  !refference (first order)     
      double precision ecorr(nstat,nstat) !correction 2. to wanted order
      double precision e(nstat,nstat)     !full diabatic matrix

      double precision q(qn)   !coordinates
      integer npar             !number of parameters
      integer pst(2,np)        !pointers for parameters
      double precision p(npar) !parameters

      double precision u(nstat,nstat)    !ci-vectors
      double precision uref(nstat,nstat) !ci-vectors
      double precision v(nstat)          !eigen values
      double precision dv(nstat)          !eigen values
      double precision vref(nstat)       !eigen values

      double precision damp !for fixing black holes

c     Parameters for damping function
      double precision minimum
      double precision width
      parameter(minimum = -0.2, width =0.1)

      logical skip

      skip=.false.

      n=1

c     compute reference model:
      iref=0
      call diab(eref,n,q,p,pst,iref,npar)

c     compute correction  model:
      iref=1
      call diab(ecorr,n,q,p,pst,iref,npar)

c     sum refference and correction
      if (skip) return
      do i=1,nstat
        do j=1,nstat
          e(i,j)=eref(i,j)+ecorr(i,j)
        enddo
      enddo

      call rs(nstat,nstat,eref,vref,1,uref)
      call rs(nstat,nstat,e,v,1,u)

      do i=1,3
         dv(i) = v(i) - vref(i)
      enddo

c     now compute full model with damped correction matrix:
      damp = ((tanh((dv(1)-minimum)/width) + 1.)/2.)

      do i=1,nstat
        do j=1,nstat
          e(i,j)=eref(i,j)+damp*ecorr(i,j)
        enddo
      enddo

      return

      end
