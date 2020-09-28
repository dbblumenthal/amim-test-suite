package de.mpg.mpiinf.ag1.kpm.utils;

/**
 *
 * @author nalcaraz
 */
public enum Separator {
    TAB, SPACE, COMMA;
    
    public char charValue() {
        switch (this) {
            case TAB:
                return '\t';
            case SPACE:
                return ' ';
            case COMMA:
                return ',';
        }
        return '\t';       
    }
}
