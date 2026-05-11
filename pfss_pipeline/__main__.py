import os

# Force certifi's CA bundle. This conda env's `base` activate hook exports
# SSL_CERT_FILE=/usr/share/ssl/certs/ca-bundle.crt, a 2019-era system bundle
# that lacks Let's Encrypt's ISRG Root X1 — JSOC HTTPS fails against it.
# DEM's own activate.d/ssl_cert.sh sets certifi correctly, but only fires on
# `conda activate DEM`, not on direct PATH-based python invocations.
try:
    import certifi
    os.environ["SSL_CERT_FILE"] = certifi.where()
    os.environ["REQUESTS_CA_BUNDLE"] = certifi.where()
except ImportError:
    pass

from pfss_pipeline.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
