import java.math.BigInteger;
import java.util.*;
import java.io.BufferedReader;
import java.io.InputStreamReader;

public class LagrangePolynomial {

    // -------- Fraction Class --------
    static class Fraction {
        BigInteger num;
        BigInteger den;

        Fraction(BigInteger n, BigInteger d) {
            if (d.signum() == 0)
                throw new ArithmeticException("Denominator zero");

            if (d.signum() < 0) {
                n = n.negate();
                d = d.negate();
            }

            BigInteger gcd = n.gcd(d);
            num = n.divide(gcd);
            den = d.divide(gcd);
        }

        static Fraction ZERO = new Fraction(BigInteger.ZERO, BigInteger.ONE);
        static Fraction ONE = new Fraction(BigInteger.ONE, BigInteger.ONE);

        Fraction add(Fraction f) {
            BigInteger n = this.num.multiply(f.den).add(f.num.multiply(this.den));
            BigInteger d = this.den.multiply(f.den);
            return new Fraction(n, d);
        }

        Fraction multiply(Fraction f) {
            return new Fraction(this.num.multiply(f.num), this.den.multiply(f.den));
        }

        Fraction subtract(Fraction f) {
            return add(new Fraction(f.num.negate(), f.den));
        }

        Fraction divide(Fraction f) {
            return new Fraction(this.num.multiply(f.den), this.den.multiply(f.num));
        }

        public String toString() {
            if (den.equals(BigInteger.ONE))
                return num.toString();
            return num + "/" + den;
        }
    }

    // Convert base string to BigInteger
    static BigInteger convertBase(String value, int base) {
        return new BigInteger(value, base);
    }

    // Simple JSON parsing without external library
    static class SimpleJSON {
        static Map<String, Object> parse(String json) {
            Map<String, Object> result = new LinkedHashMap<>();
            json = json.trim();
            if (!json.startsWith("{") || !json.endsWith("}")) {
                throw new RuntimeException("Invalid JSON");
            }
            
            int braceCount = 0;
            int start = 1;
            boolean inString = false;
            boolean escape = false;
            
            for (int i = 1; i < json.length() - 1; i++) {
                char c = json.charAt(i);
                
                if (escape) {
                    escape = false;
                    continue;
                }
                
                if (c == '\\') {
                    escape = true;
                    continue;
                }
                
                if (c == '"') {
                    inString = !inString;
                    continue;
                }
                
                if (inString) continue;
                
                if (c == '{') braceCount++;
                if (c == '}') braceCount--;
                
                if (c == ',' && braceCount == 0) {
                    String pair = json.substring(start, i).trim();
                    addPair(result, pair);
                    start = i + 1;
                }
            }
            
            String lastPair = json.substring(start, json.length() - 1).trim();
            if (!lastPair.isEmpty()) {
                addPair(result, lastPair);
            }
            
            return result;
        }
        
        private static void addPair(Map<String, Object> map, String pair) {
            int colonIndex = pair.indexOf(':');
            if (colonIndex == -1) return;
            
            String key = pair.substring(0, colonIndex).trim();
            String value = pair.substring(colonIndex + 1).trim();
            
            // Remove quotes from key
            if (key.startsWith("\"") && key.endsWith("\"")) {
                key = key.substring(1, key.length() - 1);
            }
            
            // Parse value
            Object parsedValue;
            if (value.startsWith("{")) {
                parsedValue = parse(value);
            } else if (value.startsWith("\"")) {
                parsedValue = value.substring(1, value.length() - 1);
            } else if (value.equals("true")) {
                parsedValue = Boolean.TRUE;
            } else if (value.equals("false")) {
                parsedValue = Boolean.FALSE;
            } else if (value.equals("null")) {
                parsedValue = null;
            } else {
                try {
                    parsedValue = Integer.parseInt(value);
                } catch (NumberFormatException e) {
                    try {
                        parsedValue = Long.parseLong(value);
                    } catch (NumberFormatException e2) {
                        parsedValue = value;
                    }
                }
            }
            
            map.put(key, parsedValue);
        }
        
        @SuppressWarnings("unchecked")
        static Map<String, Object> getObject(Map<String, Object> map, String key) {
            Object obj = map.get(key);
            if (obj instanceof Map) {
                return (Map<String, Object>) obj;
            }
            throw new RuntimeException("Key " + key + " is not an object");
        }
        
        static int getInt(Map<String, Object> map, String key) {
            Object obj = map.get(key);
            if (obj instanceof Number) {
                return ((Number) obj).intValue();
            }
            return Integer.parseInt(obj.toString());
        }
        
        static String getString(Map<String, Object> map, String key) {
            Object obj = map.get(key);
            return obj.toString();
        }
    }

    public static void main(String[] args) throws Exception {

        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        StringBuilder jsonInput = new StringBuilder();
        String line;

        while ((line = br.readLine()) != null) {
            jsonInput.append(line);
        }

        Map<String, Object> json = SimpleJSON.parse(jsonInput.toString());

        Map<String, Object> keys = SimpleJSON.getObject(json, "keys");
        int k = SimpleJSON.getInt(keys, "k");

        // Collect all roots and sort by x value
        List<Map.Entry<BigInteger, BigInteger>> roots = new ArrayList<>();

        for (String key : json.keySet()) {
            if (key.equals("keys")) continue;

            BigInteger x = new BigInteger(key);
            Map<String, Object> root = SimpleJSON.getObject(json, key);

            int base = Integer.parseInt(SimpleJSON.getString(root, "base"));
            String value = SimpleJSON.getString(root, "value");

            BigInteger y = convertBase(value, base);

            roots.add(Map.entry(x, y));
        }

        // Sort by x value
        roots.sort(Map.Entry.comparingByKey());

        // Use first k points
        int pointsToUse = Math.min(k, roots.size());

        // Lagrange Interpolation to find f(0)
        // f(0) = Σ yi * Li(0)
        // Li(0) = Π (0 - xj) / (xi - xj) for j ≠ i

        Fraction constantTerm = Fraction.ZERO;

        for (int i = 0; i < pointsToUse; i++) {
            BigInteger xi = roots.get(i).getKey();
            BigInteger yi = roots.get(i).getValue();

            // Calculate Lagrange basis polynomial L_i(0)
            Fraction li = Fraction.ONE;

            for (int j = 0; j < pointsToUse; j++) {
                if (i == j) continue;

                BigInteger xj = roots.get(j).getKey();

                // L_i(0) = L_i(0) * (0 - xj) / (xi - xj)
                //        = L_i(0) * (-xj) / (xi - xj)
                Fraction numerator = new Fraction(xj.negate(), BigInteger.ONE);
                Fraction denominator = new Fraction(xi.subtract(xj), BigInteger.ONE);

                li = li.multiply(numerator).divide(denominator);
            }

            // Add yi * Li(0) to the result
            Fraction yiFraction = new Fraction(yi, BigInteger.ONE);
            constantTerm = constantTerm.add(yiFraction.multiply(li));
        }

        System.out.println(constantTerm);
    }
}
