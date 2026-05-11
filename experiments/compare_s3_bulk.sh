set -euo pipefail
cd /Users/spuduch/Research/leiomyoma_proteomics

# Set creds in env (avoid putting secrets in command history)
export AWS_ACCESS_KEY_ID='SE0AV795UKCQ338YKWP4'
export AWS_SECRET_ACCESS_KEY='/mkkvYtFJkO+NAhxcm3OhNKAdvwQivhbdQRLeJ/c'
export AWS_DEFAULT_REGION='us-east-1'

ENDPOINT='https://s3-ext.decode.is:10443'
BUCKET='largescaleplasma-2023'
N=40
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# Same sample for both methods
awk 'index($0,"file=") && $0 ~ /\.txt\.gz([[:space:]]*)$/ {print}' data/raw/deCODE/bulk_urls.txt \
| awk 'BEGIN{srand()} {print rand() "\t" $0}' \
| sort -k1,1n \
| awk -F '\t' -v n="$N" 'NR <= n {print $2}' > "$TMP/http_urls.txt"

echo "Sample URL count: $(wc -l < "$TMP/http_urls.txt")"
if [ ! -s "$TMP/http_urls.txt" ]; then
  echo "No .txt.gz URLs matched in data/raw/deCODE/bulk_urls.txt"
  exit 1
fi

sed -E 's|.*file=([^&]+).*|\1|' "$TMP/http_urls.txt" > "$TMP/keys.txt"
echo "Sample key count: $(wc -l < "$TMP/keys.txt")"

echo "Identity check on 5 files (HTTP vs S3 SHA256)"
awk 'NR <= 5' "$TMP/http_urls.txt" > "$TMP/http_urls_5.txt"
awk 'NR <= 5' "$TMP/keys.txt" > "$TMP/keys_5.txt"
paste "$TMP/http_urls_5.txt" "$TMP/keys_5.txt" | while IFS=$'\t' read -r url key; do
  h1=$(curl -fsSL "$url" | shasum -a 256 | awk '{print $1}')
  h2=$(aws s3 cp "s3://$BUCKET/$key" - --endpoint-url "$ENDPOINT" --no-progress 2>/dev/null | shasum -a 256 | awk '{print $1}')
  [[ "$h1" == "$h2" ]] && echo "OK $key" || echo "MISMATCH $key"
done

echo
echo "HTTP benchmark (workers 4,12)"
for w in 4 12; do
  echo "workers=$w"
  /usr/bin/time -p sh -c \
    "cat '$TMP/http_urls.txt' | xargs -n1 -P $w -I{} curl -fsSL '{}' -o /dev/null"
done

echo
echo "S3 benchmark (workers 4,12)"
for w in 4 12; do
  echo "workers=$w"
  /usr/bin/time -p sh -c \
    "cat '$TMP/keys.txt' | xargs -n1 -P $w -I{} aws s3 cp 's3://$BUCKET/{}' - --endpoint-url '$ENDPOINT' --no-progress >/dev/null"
done
