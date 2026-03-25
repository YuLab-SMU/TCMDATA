#!/usr/bin/env Rscript
# Local integration test for download_ligand_structure,
# download_receptor_structure, and convert_structure.

devtools::load_all(".")

tmpdir <- file.path(tempdir(), "tcmdata_test")
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
cat("Test output dir:", tmpdir, "\n\n")

passed <- 0
failed <- 0

assert <- function(desc, expr) {
  ok <- tryCatch(expr, error = function(e) { message("  ERROR: ", e$message); FALSE })
  if (isTRUE(ok)) {
    cat("[PASS]", desc, "\n")
    passed <<- passed + 1
  } else {
    cat("[FAIL]", desc, "\n")
    failed <<- failed + 1
  }
}

# ── 0. Print obabel location ──────────────────────────────────────
cat("== OpenBabel location ==\n")
obabel_path <- Sys.which("obabel")
cat("obabel:", obabel_path, "\n")
system2(obabel_path, "-V", stdout = TRUE, stderr = TRUE) |> cat(sep = "\n")
cat("\n")

# ── 1. download_ligand_structure ──────────────────────────────────
cat("== Test: download_ligand_structure ==\n")

sdf_file <- download_ligand_structure(2244, destdir = tmpdir, overwrite = TRUE)
assert("returns a file path", is.character(sdf_file))
assert("SDF file exists", file.exists(sdf_file))
assert("SDF contains $$$$ marker",
       any(grepl("\\$\\$\\$\\$", readLines(sdf_file))))
cat("\n")

# ── 2. download_receptor_structure ────────────────────────────────
cat("== Test: download_receptor_structure ==\n")

pdb_file <- download_receptor_structure("1crn", destdir = tmpdir, overwrite = TRUE)
assert("returns a file path", is.character(pdb_file))
assert("PDB file exists", file.exists(pdb_file))
assert("PDB contains ATOM records",
       any(grepl("^ATOM", readLines(pdb_file))))

cif_file <- download_receptor_structure("1crn", format = "cif",
                                        destdir = tmpdir, overwrite = TRUE)
assert("CIF file exists", file.exists(cif_file))
assert("CIF contains data_ header",
       any(grepl("^data_", readLines(cif_file))))
cat("\n")

# ── 3. convert_structure ──────────────────────────────────────────
cat("== Test: convert_structure ==\n")

# helper: build output path in tmpdir
out <- function(name) file.path(tmpdir, name)

# 3a. SDF -> PDB
pdb_conv <- convert_structure(sdf_file, output_file = out("aspirin.pdb"),
                              output_type = "pdb", overwrite = TRUE)
assert("SDF->PDB: output exists", file.exists(pdb_conv))
pdb_conv_lines <- readLines(pdb_conv)
assert("SDF->PDB: contains ATOM/HETATM",
       any(grepl("^(ATOM|HETATM)", pdb_conv_lines)))

# 3b. SDF -> MOL2
mol2_conv <- convert_structure(sdf_file, output_file = out("aspirin.mol2"),
                               output_type = "mol2", overwrite = TRUE)
assert("SDF->MOL2: output exists", file.exists(mol2_conv))
assert("SDF->MOL2: contains @<TRIPOS>MOLECULE",
       any(grepl("@<TRIPOS>MOLECULE", readLines(mol2_conv))))

# 3c. PDB -> XYZ
xyz_conv <- convert_structure(pdb_conv, output_file = out("aspirin.xyz"),
                              output_type = "xyz", overwrite = TRUE)
assert("PDB->XYZ: output exists", file.exists(xyz_conv))
xyz_lines <- readLines(xyz_conv)
assert("PDB->XYZ: first line is atom count (integer)",
       grepl("^\\s*\\d+\\s*$", xyz_lines[1]))

# 3d. SDF -> SMI (SMILES)
smi_conv <- convert_structure(sdf_file, output_file = out("aspirin.smi"),
                              output_type = "smi", overwrite = TRUE)
assert("SDF->SMI: output exists", file.exists(smi_conv))
assert("SDF->SMI: non-empty output",
       length(readLines(smi_conv)) > 0 && nchar(readLines(smi_conv)[1]) > 0)

# 3e. Explicit input_type
pdb_conv2 <- convert_structure(sdf_file, input_type = "sdf", output_type = "pdb",
                               output_file = out("explicit_test.pdb"),
                               overwrite = TRUE)
assert("Explicit input_type: output exists", file.exists(pdb_conv2))

# 3f. add_hydrogens produces >= atoms
pdb_h <- convert_structure(sdf_file, output_type = "pdb",
                           add_hydrogens = TRUE,
                           output_file = out("with_h.pdb"),
                           overwrite = TRUE)
assert("add_hydrogens: output exists", file.exists(pdb_h))
h_lines <- readLines(pdb_h)
assert("add_hydrogens: atom count >= without hydrogens",
       sum(grepl("^(ATOM|HETATM)", h_lines)) >=
         sum(grepl("^(ATOM|HETATM)", pdb_conv_lines)))

cat("\n== Results ==\n")
cat(sprintf("Passed: %d / %d\n", passed, passed + failed))
if (failed > 0) {
  cat(sprintf("FAILED: %d\n", failed))
  quit(status = 1)
} else {
  cat("All tests passed!\n")
}
