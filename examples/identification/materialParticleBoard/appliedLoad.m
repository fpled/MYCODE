function F = appliedLoad(sample)
% function F = appliedLoad(sample)
% input :
% sample = 'B1',...,'B27'   Elastic beam sample 
%          'S1',...,'S16'   Screw junction sample    
%          'D1',...,'D8'    Dowel junction sample
% output :
% F = applied load [N]

scal_DELTALAB = 0.7;
scal_INSTRON = 0.55;

switch sample
    
    case 'B1'
        F = [48 113 157 217 253 313 351 413 453 502]*scal_DELTALAB;
    case 'B2'
        F = [10 70 111 169 209 269 309 347 399]*scal_DELTALAB;
    case 'B3'
        F = [18 54 115 156 221 258 314 360 413 452]*scal_DELTALAB;
    case 'B4'
        F = [13 68  98 168 209 262 301 360 399 458]*scal_DELTALAB;
    case 'B5'
        F = [19 50 117 162 222 259 316 353 413 453]*scal_DELTALAB;
    case 'B6'
        F = [17 63 109 158 214 259 318 348 414 450]*scal_DELTALAB;
    case 'B7'
        F = [15 55 122 162 216 256 315 358 415 450]*scal_DELTALAB;
    case 'B8'
        F = [20 60 119 159 219 255 311 356 414 448]*scal_DELTALAB;
    case 'B9'
        F = [17 59 117 161 215 258 311 361 411 445]*scal_DELTALAB;
    case 'B10'
        F = [19 58 120 159 223 257 320 357 414 450]*scal_DELTALAB;
    case 'B11'
        F = [ 6 63 108 160 211 269 310 366 405 461 500]*scal_DELTALAB;
    case 'B12'
        F = [19 59 115 159 219 260 315 350]*scal_DELTALAB;
    case 'B13'
        F = [ 9 69  98 169 208 269 318 366 409 463]*scal_DELTALAB;
    case 'B14'
        F = [22 54 120 161 218 255 307 359 395 455]*scal_DELTALAB;
    case 'B15'
        F = [19 66 108 165 207 264 304 358 392]*scal_DELTALAB;
    case 'B16'
        F = [ 8 73 107 168 204 273 310 365 403 470]*scal_DELTALAB;
    case 'B17'
        F = [49 99 148 198 245 295 345 395 442 490 543 585]*scal_INSTRON;
    case 'B18'
        F = [48 98 147 197 247 296 346 393 442 490 538 586]*scal_INSTRON;
    case 'B19'
        F = [48 98 147 197 247 297 347 395 445 490 540]*scal_INSTRON;
    case 'B20'
        F = [48 98 148 197 247 296 344 393 443 492 540]*scal_INSTRON;
    case 'B21'
        F = [48 98 148 196 246 295 345 390 440 490 535]*scal_INSTRON;
    case 'B22'
        F = [46 95 141 193 242 291 338 386 430 477 430]*scal_INSTRON;
    case 'B23'
        F = [45 96 146 195 244 292 340 390 438 482 515 574]*scal_INSTRON;
    case 'B24'
        F = [47 97 146 195 245 295 340 390 440 490 530]*scal_INSTRON;
    case 'B25'
        F = [47 97 147 196 245 293 340 390 438 488 525]*scal_INSTRON;
    case 'B26'
        F = [17 68 118 159 208 257 313 355 410 447]*scal_DELTALAB;
    case 'B27'
        F = [19 57 114 159 218 257 316 355 393 366]*scal_DELTALAB;
        
        
    case 'S1'
        F = [24 65 114 163 200 252 291];
    case 'S2'
        F = [14 54 113 159 199 249 299];
    case 'S3'
        F = [15 65 111 153 204 243 295 341 389];
    case 'S4'
        F = [23 34 91 141 180 229 290];
    case 'S5'
        F = [41 89 128 182 250 291 330];
    case 'S6'
        F = [33 81 123 184 251 294 334 380];
    case 'S7'
        F = [26 87 123 175 241 287 329];
    case 'S8'
        F = [42 92 140 189 254 298 341 384];
    case 'S9'
        F = [24 74 121 169 231 256];
    case 'S10'
        F = [26 68 110 156 218];
    case 'S11'
        F = [33 81 133 177 202 245 289];
    case 'S12'
        F = [36 79 126 149 187 241 292 329];
    case 'S13'
        F = [30 80 126 177 244 287 330];
    case 'S14'
        F = [41 82 135 189 246 288];
    case 'S15'
        F = [27 75 123 167 238 266];
    case 'S16'
        F = [25 74 118 173 218];
        
        
    case 'D1'
        % F = [12 24 98 146];
        F = [24 49 98 146];
    case 'D2'
        % F = [6 63 110 152 197];
        F = [12 24 110 152 98];
    case 'D3'
        % F = [12 38 48 61];
        F = [12 37 49 61];
    case 'D4'
        % F = [12 36];
        F = [12 37];
    case 'D5'
        % F = [12 24 36];
        F = [12 24 37];
    case 'D6'
        % F = [12 35 48];
        F = [12 37 49];
    case 'D7'
        % F = [12 38 61];
        F = [12 37 61];
    case 'D8'
        % F = [12 36 48];
        F = [12 37 49];
        
    otherwise
        error(['Wrong sample' sample])
        
end

