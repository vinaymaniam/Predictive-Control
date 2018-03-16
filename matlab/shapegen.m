function [sp, tp, c] = shapegen(i)
    if i > 7
        i = 5;
    end
    switch i
        case 1
            sp = [0.05, 0.4];
            tp = [0.4, 0.05];
            mp = [0.4, 0.4];
            c = round(shape(sp, tp, 0.1, mp, 0),4);
        case 2
            tp = [0.05, 0.4];
            sp = [0.4, 0.05];
            mp = [0.05, 0.05];
            c = round(shape(sp, tp, 0.1, mp, 0),4);
        case 3
            sp = [0.4, 0.4];
            tp = [0.05, 0.05];
            mp = [0.4, 0.05];
            c = round(shape(sp, tp, 0.1, mp, 0),4);
        case 4
            sp = [0.05, 0.05];
            tp = [0.4, 0.4];
            mp = [0.05, 0.4];
            c = round(shape(sp, tp, 0.1, mp, 0),4);
        case 5        
            sp = [0.05, 0.05];
            tp = [0.45, 0.05];
            mp = [0.25, 0.25];
            c = round(shape(sp, tp, 0.1, mp, 0),4);
        case 6
            sp = [0.05, 0.45];
            tp = [0.05, 0.05];
            mp = [0.25, 0.25];
            c = round(shape(sp, tp, 0.1, mp, 0),4);
        case 7
            sp = [0.45, 0.3];
            tp = [0.05, 0.3];
            mp = [0.25, 0.05];
            c = round(shape(sp, tp, 0.1, mp, 0),4);        
    end
end