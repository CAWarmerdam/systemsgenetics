package nl.systemsgenetics.pgsbasedmixupmapper;

public class PGSBasedMixupMapperException extends Exception {
    public PGSBasedMixupMapperException(String s, Exception e) {
        super(s, e);
    }

    public PGSBasedMixupMapperException(String s) {
        super(s);
    }
}