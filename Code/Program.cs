using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

// Run this:
//
// dotnet new console -n CollatzCert
// cd CollatzCert
// # Replace Program.cs with this file
// dotnet run --file phi_k13_conservative_certificate.csv --top 10 --report report.csv --startN 513
// # Stronger policy (force conservative exceptional valuation):
// dotnet run --file phi_k13_conservative_certificate.csv --force-exceptional

class Program
{
    // === Verification policy flag (now settable via CLI) =================
    // Use --force-exceptional to *enforce* the conservative valuation on
    // the exceptional residue r* (the unique odd r with (3*r+1) % 2^k == 0):
    // namely v2 := k and the canonical target F_k(r*) = odd(3r*+1) mod 2^k.
    // This makes the check strictly harder and documents the intended policy.
    static bool ForceConservativeExceptional = false;
    // =====================================================================

    // ===== Utils =====
    static int V2(long n) { int c = 0; while ((n & 1) == 0) { n >>= 1; c++; } return c; }
    static long OddPart(long n) { while ((n & 1) == 0) n >>= 1; return n; }
    static int InferKFromResidues(IEnumerable<int> residues)
    {
        int rmax = residues.Max();
        int k = (int)Math.Round(Math.Log(rmax + 1, 2.0));
        if ((1 << k) != rmax + 1)
        {
            int n = residues.Count();
            k = (int)Math.Round(Math.Log(n * 2.0, 2.0));
        }
        return k;
    }

    static List<string[]> ReadCsv(string path)
    {
        var lines = File.ReadAllLines(path);
        var rows = new List<string[]>();
        foreach (var line in lines) rows.Add(ParseCsvLine(line));
        return rows;
    }
    static string[] ParseCsvLine(string line)
    {
        var cells = new List<string>();
        bool inQuotes = false;
        var cur = new System.Text.StringBuilder();
        for (int i = 0; i < line.Length; i++)
        {
            char c = line[i];
            if (inQuotes)
            {
                if (c == '"')
                {
                    if (i + 1 < line.Length && line[i + 1] == '"') { cur.Append('"'); i++; }
                    else inQuotes = false;
                }
                else cur.Append(c);
            }
            else
            {
                if (c == ',') { cells.Add(cur.ToString()); cur.Clear(); }
                else if (c == '"') inQuotes = true;
                else cur.Append(c);
            }
        }
        cells.Add(cur.ToString());
        return cells.Select(s => s.Trim()).ToArray();
    }

    class Entry
    {
        public int r;      // odd residue mod 2^k
        public double phi; // phi(r)
        public int v2;     // v2(3r+1)
        public int fk;     // F_k(r) = odd(3r+1) mod 2^k
        public double slackFloat;
        public double slackUpper; // interval-style upper slack (must be <= 0)
    }

    static int IndexOfAny(string[] header, string[] candidates)
    {
        for (int j = 0; j < header.Length; j++)
            foreach (var c in candidates)
                if (header[j].Contains(c.ToLowerInvariant())) return j;
        return -1;
    }

    static void Main(string[] args)
    {
        // ===== Parameters =====
        string path = "phi_k13_conservative_certificate.csv";
        double delta = 0.10, b = 0.05, p = 0.05, zeta = 0.02;
        long startN = -1;
        int top = 10;
        string reportPath = ""; // if set, write report CSV

        // CLI
        for (int i = 0; i < args.Length; i++)
        {
            string a = args[i];
            if (a == "--file" && i + 1 < args.Length) path = args[++i];
            else if (a == "--delta" && i + 1 < args.Length) delta = double.Parse(args[++i], CultureInfo.InvariantCulture);
            else if (a == "--b" && i + 1 < args.Length) b = double.Parse(args[++i], CultureInfo.InvariantCulture);
            else if (a == "--p" && i + 1 < args.Length) p = double.Parse(args[++i], CultureInfo.InvariantCulture);
            else if (a == "--zeta" && i + 1 < args.Length) zeta = double.Parse(args[++i], CultureInfo.InvariantCulture);
            else if (a == "--startN" && i + 1 < args.Length) startN = long.Parse(args[++i], CultureInfo.InvariantCulture);
            else if (a == "--top" && i + 1 < args.Length) top = int.Parse(args[++i], CultureInfo.InvariantCulture);
            else if (a == "--report" && i + 1 < args.Length) reportPath = args[++i];
            else if (a == "--force-exceptional") ForceConservativeExceptional = true;
        }

        if (!File.Exists(path))
        {
            Console.WriteLine($"[ERROR] CSV not found: {path}");
            Environment.Exit(10);
        }

        // ===== Load CSV =====
        var rows = ReadCsv(path);
        if (rows.Count == 0) { Console.WriteLine("[ERROR] Empty CSV"); Environment.Exit(11); }

        var header = rows[0].Select(s => s.ToLowerInvariant()).ToArray();
        int idxR   = IndexOfAny(header, new[] { "r", "residue", "odd_residue", "odd_residue_mod_2^k" });
        int idxPhi = IndexOfAny(header, new[] { "phi", "ϕ", "varphi" });
        int idxV2  = IndexOfAny(header, new[] { "v2", "nu2", "v2(3r+1)", "valuation", "nu_2", "nu-2" });
        int idxFk  = IndexOfAny(header, new[] { "fk", "f_k", "next", "target" });
        if (idxR < 0 || idxPhi < 0)
        {
            Console.WriteLine("[ERROR] Missing residue/phi columns.");
            Console.WriteLine("Header: " + string.Join(" | ", rows[0]));
            Environment.Exit(12);
        }

        var data = new List<Entry>();
        for (int i = 1; i < rows.Count; i++)
        {
            var row = rows[i];
            if (row.Length == 0 || string.IsNullOrWhiteSpace(row[0])) continue;
            int r = int.Parse(row[idxR], CultureInfo.InvariantCulture);
            double phi = double.Parse(row[idxPhi], CultureInfo.InvariantCulture);
            int v2 = -1, fk = -1;
            if (idxV2 >= 0 && idxV2 < row.Length && !string.IsNullOrWhiteSpace(row[idxV2]))
                v2 = int.Parse(row[idxV2], CultureInfo.InvariantCulture);
            if (idxFk >= 0 && idxFk < row.Length && !string.IsNullOrWhiteSpace(row[idxFk]))
                fk = int.Parse(row[idxFk], CultureInfo.InvariantCulture);
            data.Add(new Entry { r = r, phi = phi, v2 = v2, fk = fk });
        }

        int k = InferKFromResidues(data.Select(e => e.r));
        int M = 1 << k;

        foreach (var e in data)
        {
            // Compute missing values from CSV if needed
            if (e.v2 < 0) e.v2 = V2(3L * e.r + 1);
            if (e.fk < 0)
            {
                long op = OddPart(3L * e.r + 1);
                e.fk = (int)(op % M);
                if ((e.fk & 1) == 0) e.fk = (e.fk + 1) % M; // safety
            }

            // Apply *strong* policy if requested: regardless of CSV contents
            if (ForceConservativeExceptional)
            {
                bool isExceptional = ((3L * e.r + 1) % M) == 0;
                if (isExceptional)
                {
                    e.v2 = k; // minimal permissible valuation
                    long opCanon = OddPart(3L * e.r + 1);
                    int fkCanon = (int)(opCanon % M);
                    if ((fkCanon & 1) == 0) fkCanon = (fkCanon + 1) % M;
                    e.fk = fkCanon; // canonical target
                }
            }
        }

        var phiMap = data.ToDictionary(e => e.r, e => e.phi);

        // ===== Floating-point check =====
        double ln2 = Math.Log(2.0), ln3 = Math.Log(3.0);
        double rho = b * ln2 - Math.Log(1.0 - p) + zeta;
        double maxSlack = double.NegativeInfinity, minSlack = double.PositiveInfinity;
        int violFloat = 0;
        foreach (var e in data)
        {
            double phiFk = phiMap[e.fk];
            e.slackFloat = (ln3 - e.v2 * ln2) + rho + phiFk + delta - e.phi;
            if (e.slackFloat > 0) violFloat++;
            if (e.slackFloat > maxSlack) maxSlack = e.slackFloat;
            if (e.slackFloat < minSlack) minSlack = e.slackFloat;
        }

        // ===== Interval-style conservative check =====
        double L2 = 0.6931471805599452, U2 = 0.6931471805599454;
        double U3 = 1.0986122886681099;
        double kappaUpper = KappaUpperTaylor(p, 6);           // tightened
        double rhoUpper   = b * U2 + kappaUpper + zeta;       // conservative rho

        double maxSlackUp = double.NegativeInfinity, minSlackUp = double.PositiveInfinity;
        int violUp = 0;
        foreach (var e in data)
        {
            double phiFk = phiMap[e.fk];
            double lhsUpper = (U3 - e.v2 * L2) + rhoUpper + phiFk + delta;
            e.slackUpper = lhsUpper - e.phi; // must be <= 0
            if (e.slackUpper > 0) violUp++;
            if (e.slackUpper > maxSlackUp) maxSlackUp = e.slackUpper;
            if (e.slackUpper < minSlackUp) minSlackUp = e.slackUpper;
        }

        // ===== Print summary header =====
        Console.WriteLine("=== CERTIFICATE CHECK (k = {0}, |S| = {1}, M = 2^{0} = {2}) ===", k, data.Count, M);
        Console.WriteLine(ForceConservativeExceptional
            ? "Policy: --force-exceptional ENABLED (conservative exceptional valuation)"
            : "Policy: default (no forced exceptional valuation)");
        Console.WriteLine($"rho = {rho:R}, delta = {delta:R}, rho+delta = {rho + delta:R}, ln(4/3) = {Math.Log(4.0/3.0):R}");
        double margin = Math.Log(4.0 / 3.0) - (rho + delta);
        Console.WriteLine($"margin = ln(4/3) - (rho+delta) = {margin:R}");
        double Nstar = (margin > 0) ? Math.Ceiling(1.0 / (3.0 * (Math.Exp(margin) - 1.0))) : double.NaN;
        Console.WriteLine($"N* threshold (epsilon <= margin): N* = {Nstar}");

        Console.WriteLine();
        Console.WriteLine("— Floating-point check:");
        Console.WriteLine($"   max slack = {maxSlack:R}, min slack = {minSlack:R}, violations = {violFloat}");
        Console.WriteLine("— Interval-style conservative check:");
        Console.WriteLine($"   rho_upper = {rhoUpper:R} (kappaUpper={kappaUpper:R})");
        Console.WriteLine($"   max slack_upper = {maxSlackUp:R}, min slack_upper = {minSlackUp:R}, violations = {violUp}");
        Console.WriteLine();

        // ===== Optional: Top-N hardest constraints =====
        if (top > 0)
        {
            Console.WriteLine($"Top {top} (least negative) slack_upper:");
            foreach (var e in data.OrderBy(x => x.slackUpper).Take(top))
                Console.WriteLine($"  r={e.r}, fk={e.fk}, v2={e.v2}, phi(r)={e.phi:R}, slack_up={e.slackUpper:R}");
            Console.WriteLine();
        }

        // ===== Optional: write report =====
        if (!string.IsNullOrWhiteSpace(reportPath))
        {
            using var sw = new StreamWriter(reportPath);
            sw.WriteLine("r,fk,v2,phi,slack_float,slack_upper");
            foreach (var e in data.OrderBy(x => x.slackUpper))
                sw.WriteLine($"{e.r},{e.fk},{e.v2},{e.phi.ToString("R",CultureInfo.InvariantCulture)},{e.slackFloat.ToString("R",CultureInfo.InvariantCulture)},{e.slackUpper.ToString("R",CultureInfo.InvariantCulture)}");
            Console.WriteLine($"[OK] Report written: {reportPath}");
            Console.WriteLine();
        }

        // ===== Optional: simulate trajectory =====
        if (startN > 0)
        {
            if ((startN & 1) == 0) while ((startN & 1) == 0) startN >>= 1;
            Console.WriteLine($"=== Accelerated odd subsequence from N0 = {startN} ===");
            long N = startN; int step = 0;
            Console.WriteLine("t\tN_t\tv2(3N+1)\tr_t\tphi(r_t)\tlnN_t\tV_t\t→\tN_{t+1}\tr_{t+1}\tphi(r_{t+1})\tΔV\tbound");
            while (N != 1 && step < 20000)
            {
                int r = (int)(N % M); if ((r & 1) == 0) r = (r + 1) % M;
                double phi_r = phiMap[r], lnN = Math.Log(N), V = lnN + phi_r;
                long threeN1 = 3L * N + 1; int v = V2(threeN1); long Nn = OddPart(threeN1);
                int rn = (int)(Nn % M); if ((rn & 1) == 0) rn = (rn + 1) % M;
                double phi_rn = phiMap[rn], lnNn = Math.Log(Nn), Vn = lnNn + phi_rn;
                double dV = Vn - V, eps = Math.Log(1.0 + 1.0 / (3.0 * N));
                double bound = -delta + (eps - (b * ln2 - Math.Log(1.0 - p) + zeta));
                Console.WriteLine($"{step}\t{N}\t{v}\t{r}\t{phi_r:R}\t{lnN:R}\t{V:R}\t→\t{Nn}\t{rn}\t{phi_rn:R}\t{dV:R}\t{bound:R}");
                N = Nn; step++;
            }
            Console.WriteLine($"Total steps: {step}");
            Console.WriteLine();
        }

        // ===== FINAL SUMMARY & EXIT CODE =====
        bool okFloat = (violFloat == 0);
        bool okUpper = (violUp == 0);
        Console.WriteLine("=== SUMMARY ===");
        if (okFloat && okUpper)
        {
            Console.WriteLine("PASS: CERTIFICATE VERIFIED: all inequalities hold (float & interval).");
            Environment.Exit(0);
        }
        else if (!okFloat && !okUpper)
        {
            Console.WriteLine("FAIL: CERTIFICATE FAILED: violations in both float and interval checks.");
            Environment.Exit(2);
        }
        else if (!okUpper)
        {
            Console.WriteLine("FAIL: CERTIFICATE FAILED: interval-style conservative check has violations.");
            Environment.Exit(3);
        }
        else
        {
            Console.WriteLine("⚠️  Float check fails but interval-style passes (unlikely). Recheck inputs.");
            Environment.Exit(4);
        }
    }

    // Tight upper bound for kappa = -ln(1-p) via truncated series + rigorous tail
    static double KappaUpperTaylor(double p, int N = 6)
    {
        double s = 0.0, pow = p;
        for (int n = 1; n <= N; n++)
        {
            s += pow / n;
            pow *= p;
        }
        // Tail: sum_{n=N+1..∞} p^n/n <= p^{N+1} / ((N+1)*(1-p))
        double tail = pow / ((N + 1) * (1.0 - p));
        return s + tail;
    }
}
