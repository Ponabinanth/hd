/**
 * Shamir's Secret Sharing - Lagrange Interpolation
 *
 * Problem:
 * - Given n points (x, y) of a polynomial of degree k-1
 * - The y-values are encoded in different bases
 * - Some points may be "wrong" (don't lie on the true polynomial)
 * - Find the constant term f(0) using Lagrange Interpolation
 *   with majority-vote across all C(n,k) combinations
 */

const fs = require('fs');
const path = require('path');

// ─────────────────────────────────────────────
// STEP 1: Decode a value string from any base to BigInt
// (Using BigInt to handle very large numbers with exact precision)
// ─────────────────────────────────────────────
function decodeValue(base, valueStr) {
    const bigBase = BigInt(base);
    let result = 0n;
    for (const char of valueStr.toLowerCase()) {
        const digit = BigInt(parseInt(char, parseInt(base)));
        result = result * bigBase + digit;
    }
    return result;
}

// ─────────────────────────────────────────────
// STEP 2: Lagrange Interpolation to find f(0)
//
// Formula: f(0) = Σ [ y_i * Π_{j≠i} (0 - x_j) / (x_i - x_j) ]
//
// Uses BigInt + Rational arithmetic for exact (lossless) computation
// ─────────────────────────────────────────────
function lagrangeAtZero(points) {
    function gcd(a, b) {
        a = a < 0n ? -a : a;
        b = b < 0n ? -b : b;
        while (b !== 0n) { [a, b] = [b, a % b]; }
        return a;
    }
    function simplify(n, d) {
        if (d < 0n) { n = -n; d = -d; }
        const g = gcd(n < 0n ? -n : n, d);
        return [n / g, d / g];
    }
    function addFrac([n1, d1], [n2, d2]) { return simplify(n1 * d2 + n2 * d1, d1 * d2); }
    function mulFrac([n1, d1], [n2, d2]) { return simplify(n1 * n2, d1 * d2); }

    let result = [0n, 1n];
    for (let i = 0; i < points.length; i++) {
        const { x: xi, y: yi } = points[i];
        let num = 1n, den = 1n;
        for (let j = 0; j < points.length; j++) {
            if (i !== j) {
                num *= (0n - points[j].x);
                den *= (xi - points[j].x);
            }
        }
        result = addFrac(result, mulFrac([yi, 1n], simplify(num, den)));
    }
    return result[0] / result[1];
}

// ─────────────────────────────────────────────
// STEP 3: Generate all C(n, k) combinations
// ─────────────────────────────────────────────
function combinations(arr, r) {
    if (r === 0) return [[]];
    if (arr.length < r) return [];
    const [first, ...rest] = arr;
    return [
        ...combinations(rest, r - 1).map(c => [first, ...c]),
        ...combinations(rest, r)
    ];
}

// ─────────────────────────────────────────────
// STEP 4: Find the secret using majority vote
// across all C(n,k) combinations — handles bad/corrupt points
// ─────────────────────────────────────────────
function findSecret(points, k) {
    const allCombos = combinations(points, k);
    const freq = new Map();

    for (const combo of allCombos) {
        const secret = lagrangeAtZero(combo).toString();
        freq.set(secret, (freq.get(secret) || 0) + 1);
    }

    // Return the most frequent answer
    let best = null, bestCount = 0;
    for (const [val, count] of freq.entries()) {
        if (count > bestCount) {
            bestCount = count;
            best = val;
        }
    }

    return { secret: BigInt(best), count: bestCount, total: allCombos.length };
}

// ─────────────────────────────────────────────
// STEP 5: Main function — process a JSON test case
// ─────────────────────────────────────────────
function processTestCase(filePath) {
    console.log(`\n${'═'.repeat(55)}`);
    console.log(`📂 Processing: ${path.basename(filePath)}`);
    console.log('═'.repeat(55));

    const data = JSON.parse(fs.readFileSync(filePath, 'utf-8'));
    const n = data.keys.n;
    const k = data.keys.k;

    console.log(`🔑 n = ${n} (roots provided), k = ${k} (roots needed)`);
    console.log(`📐 Polynomial degree = ${k - 1}`);
    console.log(`🔀 Total combinations C(${n},${k}) = ${combinations([...Array(n)], k).length}`);

    // Decode all points
    const points = [];
    for (const key of Object.keys(data)) {
        if (key === 'keys') continue;
        const x = BigInt(key);
        const base = parseInt(data[key].base);
        const rawValue = data[key].value;
        const y = decodeValue(base, rawValue);
        console.log(`  x=${key}: base ${base}, encoded="${rawValue}" → y=${y}`);
        points.push({ x, y });
    }

    console.log(`\n� Running majority-vote across all combinations...`);
    const { secret, count, total } = findSecret(points, k);

    console.log(`\n✅ SECRET f(0) = ${secret}`);
    console.log(`   (Found in ${count} out of ${total} combinations)`);
    console.log('═'.repeat(55));

    return secret;
}

// ─────────────────────────────────────────────
// Run both test cases
// ─────────────────────────────────────────────
const tc1 = processTestCase(path.join(__dirname, 'testcase1.json'));
const tc2 = processTestCase(path.join(__dirname, 'testcase2.json'));

console.log(`\n${'★'.repeat(55)}`);
console.log(`🎯 FINAL ANSWERS:`);
console.log(`   Test Case 1 → f(0) = ${tc1}`);
console.log(`   Test Case 2 → f(0) = ${tc2}`);
console.log('★'.repeat(55));
