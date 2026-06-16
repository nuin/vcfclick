"""The single-page web UI, embedded as a string so it ships in the wheel
with no static-file packaging or JS build step. Vanilla JS talks to the
/api endpoints in app.py."""

INDEX_HTML = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>vcfclick</title>
<style>
  :root { --fg:#1f2328; --muted:#6b7280; --pri:#b45309; --bg:#faf9f7; --line:#e7e5e4; --code:#1f2328; }
  * { box-sizing: border-box; }
  body { margin:0; font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; color:var(--fg); background:var(--bg); }
  header { display:flex; align-items:center; gap:12px; padding:10px 16px; border-bottom:1px solid var(--line); background:#fff; }
  header h1 { font-size:16px; margin:0; font-weight:700; }
  header .db { font:12px ui-monospace,Menlo,monospace; color:var(--pri); background:#fff7ed; padding:2px 8px; border-radius:4px; }
  .layout { display:flex; min-height:calc(100vh - 49px); }
  aside { width:230px; border-right:1px solid var(--line); padding:12px; background:#fff; overflow:auto; }
  aside h2 { font-size:11px; text-transform:uppercase; letter-spacing:.04em; color:var(--muted); margin:0 0 6px; }
  .tbl { margin-bottom:10px; }
  .tbl > .name { font:12px ui-monospace,monospace; font-weight:600; cursor:pointer; }
  .tbl .cols { font:11px ui-monospace,monospace; color:var(--muted); padding-left:8px; display:none; }
  .tbl.open .cols { display:block; }
  main { flex:1; padding:16px; overflow:auto; }
  nav { display:flex; gap:4px; margin-bottom:14px; border-bottom:1px solid var(--line); }
  nav button { border:0; background:none; padding:8px 14px; font-size:13px; cursor:pointer; color:var(--muted); border-bottom:2px solid transparent; }
  nav button.active { color:var(--fg); border-bottom-color:var(--pri); font-weight:600; }
  .panel { display:none; }
  .panel.active { display:block; }
  textarea, input, select { font-family:ui-monospace,Menlo,monospace; font-size:13px; border:1px solid var(--line); border-radius:6px; padding:8px; width:100%; background:#fff; }
  textarea { resize:vertical; }
  .row { display:flex; gap:10px; flex-wrap:wrap; align-items:flex-end; margin-bottom:10px; }
  .row > div { flex:1; min-width:120px; }
  label { display:block; font-size:11px; color:var(--muted); margin-bottom:3px; }
  button.run { background:var(--pri); color:#fff; border:0; border-radius:6px; padding:8px 18px; font-size:13px; font-weight:600; cursor:pointer; }
  button.run:hover { background:#92400e; }
  pre.sql { background:#1f2328; color:#e7e5e4; padding:10px 12px; border-radius:6px; overflow:auto; font-size:12px; margin:10px 0; }
  .err { color:#b91c1c; background:#fef2f2; border:1px solid #fecaca; padding:8px 12px; border-radius:6px; margin:10px 0; font-size:13px; }
  .note { color:#92400e; background:#fffbeb; border:1px solid #fde68a; padding:8px 12px; border-radius:6px; margin:10px 0; font-size:13px; }
  .results { overflow:auto; border:1px solid var(--line); border-radius:6px; background:#fff; }
  table { border-collapse:collapse; width:100%; font:12px ui-monospace,monospace; }
  th, td { text-align:left; padding:5px 9px; border-bottom:1px solid var(--line); white-space:nowrap; }
  th { background:#f5f5f4; color:var(--muted); position:sticky; top:0; }
  .meta { font-size:12px; color:var(--muted); margin:6px 2px; }
  .hint { font-size:12px; color:var(--muted); }
  .grid2 { display:grid; grid-template-columns:1fr 1fr; gap:10px; }
  .chip-set { color:var(--pri); }
</style>
</head>
<body>
<header>
  <h1>vcfclick</h1>
  <span class="db" id="dbname">—</span>
  <span class="hint" style="margin-left:auto">local web UI · localhost only</span>
</header>
<div class="layout">
  <aside>
    <h2>Schema</h2>
    <div id="schema"></div>
  </aside>
  <main>
    <nav>
      <button data-tab="query" class="active">SQL</button>
      <button data-tab="ask">Ask</button>
      <button data-tab="trio">Trio</button>
      <button data-tab="combine">Combine</button>
    </nav>

    <section class="panel active" id="panel-query">
      <textarea id="q-sql" rows="5">SELECT chrom, pos, ref, alt FROM variants LIMIT 20</textarea>
      <div style="margin:8px 0"><button class="run" onclick="runQuery()">Run</button></div>
      <div id="q-out"></div>
    </section>

    <section class="panel" id="panel-ask">
      <div class="row">
        <div style="max-width:140px">
          <label>provider</label>
          <select id="a-provider"><option value="gemini">Gemini</option><option value="anthropic">Anthropic</option></select>
        </div>
        <div><label>API key (kept in your browser)</label><input id="a-key" type="password" placeholder="paste key" /></div>
        <div style="max-width:200px"><label>model (optional)</label><input id="a-model" placeholder="default" /></div>
      </div>
      <textarea id="a-q" rows="2" placeholder="Which variants have the highest allele frequency?"></textarea>
      <div style="margin:8px 0"><button class="run" onclick="askNl()">Ask</button></div>
      <div id="a-out"></div>
    </section>

    <section class="panel" id="panel-trio">
      <div class="row">
        <div style="max-width:220px"><label>proband</label><select id="t-proband"></select></div>
        <div style="max-width:150px"><label>category</label>
          <select id="t-cat"><option>denovo</option><option>recessive</option><option>dominant</option></select></div>
        <div style="max-width:90px"><label>min GQ</label><input id="t-gq" type="number" value="20" /></div>
        <div style="max-width:90px"><label>min DP</label><input id="t-dp" type="number" value="10" /></div>
        <div style="max-width:90px"><label>max AF</label><input id="t-af" type="number" step="0.01" value="0.01" /></div>
      </div>
      <div style="margin:8px 0"><button class="run" onclick="runTrio()">Run</button></div>
      <div id="t-out"></div>
    </section>

    <section class="panel" id="panel-combine">
      <p class="hint">Combine two call sets that may share samples (GATK3 CombineVariants). Input order is priority.</p>
      <div class="grid2">
        <div><label>first.vcf</label><textarea id="c-a" rows="7">##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	101	.	A	T	.	.	.	GT	0/1	0/1
chr1	202	.	C	A	.	.	.	GT	1/1	0/1</textarea></div>
        <div><label>second.vcf</label><textarea id="c-b" rows="7">##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S2	S3
chr1	202	.	C	A	.	.	.	GT	1/1	0/1
chr1	404	.	T	G	.	.	.	GT	0/0	1/1</textarea></div>
      </div>
      <div class="row" style="margin-top:10px">
        <div style="max-width:130px"><label>--min-callsets</label><input id="c-min" type="number" value="1" min="1" max="2" /></div>
        <div style="flex:0"><button class="run" onclick="runCombine()">Combine</button></div>
      </div>
      <div id="c-out"></div>
    </section>
  </main>
</div>
<script>
const $ = (id) => document.getElementById(id);

function esc(s){ return String(s).replace(/[&<>]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;'}[c])); }

function renderTable(cols, rows){
  if(!rows || !rows.length) return '<div class="meta">0 rows</div>';
  let h = '<div class="results"><table><thead><tr>' + cols.map(c=>'<th>'+esc(c)+'</th>').join('') + '</tr></thead><tbody>';
  for(const r of rows){ h += '<tr>' + r.map(v=>'<td>'+esc(v===null?'∅':v)+'</td>').join('') + '</tr>'; }
  return h + '</tbody></table></div><div class="meta">'+rows.length+' rows</div>';
}
function showResult(el, d){
  let h = '';
  if(d.sql) h += '<pre class="sql">'+esc(d.sql)+'</pre>';
  if(d.error){ h += '<div class="err">'+esc(d.error)+'</div>'; }
  if(d.note){ h += '<div class="note">'+esc(d.note)+'</div>'; }
  if(d.columns) h += renderTable(d.columns, d.rows);
  el.innerHTML = h || '<div class="meta">no result</div>';
}
async function postJson(url, body){
  const r = await fetch(url, {method:'POST', headers:{'content-type':'application/json'}, body:JSON.stringify(body)});
  return r.json();
}

async function loadMeta(){
  const m = await (await fetch('/api/meta')).json();
  $('dbname').textContent = m.db || '(default)';
  $('schema').innerHTML = m.tables.map(t =>
    '<div class="tbl"><div class="name" onclick="this.parentNode.classList.toggle(\'open\')">'+esc(t.name)+
    '</div><div class="cols">'+t.columns.map(esc).join('<br>')+'</div></div>').join('') || '<div class="hint">empty database</div>';
  $('t-proband').innerHTML = (m.probands.length?m.probands:m.samples).map(s=>'<option>'+esc(s)+'</option>').join('');
}

async function runQuery(){
  showResult($('q-out'), await postJson('/api/query', {sql:$('q-sql').value}));
}
async function askNl(){
  const key = $('a-key').value; localStorage.setItem('vcfclick_key', key);
  localStorage.setItem('vcfclick_provider', $('a-provider').value);
  $('a-out').innerHTML = '<div class="meta">thinking…</div>';
  showResult($('a-out'), await postJson('/api/nl', {
    question:$('a-q').value, provider:$('a-provider').value, key, model:$('a-model').value }));
}
async function runTrio(){
  const p = new URLSearchParams({proband:$('t-proband').value, category:$('t-cat').value,
    min_gq:$('t-gq').value, min_dp:$('t-dp').value, max_af:$('t-af').value});
  showResult($('t-out'), await (await fetch('/api/trio?'+p)).json());
}
async function runCombine(){
  const d = await postJson('/api/combine', {first:$('c-a').value, second:$('c-b').value, min_callsets:Number($('c-min').value)});
  if(d.error){ $('c-out').innerHTML = '<div class="err">'+esc(d.error)+'</div>'; return; }
  const cols = ['pos','set=', ...d.samples];
  const rows = d.records.map(r => [r.chrom+':'+r.pos, r.set, ...d.samples.map(s=>r.cells[s]||'./.')]);
  $('c-out').innerHTML = renderTable(cols, rows);
}

document.querySelectorAll('nav button').forEach(b => b.onclick = () => {
  document.querySelectorAll('nav button').forEach(x=>x.classList.remove('active'));
  document.querySelectorAll('.panel').forEach(x=>x.classList.remove('active'));
  b.classList.add('active'); $('panel-'+b.dataset.tab).classList.add('active');
});
$('a-key').value = localStorage.getItem('vcfclick_key') || '';
if(localStorage.getItem('vcfclick_provider')) $('a-provider').value = localStorage.getItem('vcfclick_provider');
loadMeta();
</script>
</body>
</html>
"""
