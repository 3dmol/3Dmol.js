// Build first; start learning-environment/server.py --port 5072 separately.
// Override LEARNING_TEST_ORIGIN when testing another local server.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const puppeteer = require('puppeteer');
const origin = process.env.LEARNING_TEST_ORIGIN || 'http://127.0.0.1:5072';
const room = 'embed-test-' + Date.now();
(async () => {
    const browser = await puppeteer.launch({headless: true, args: ['--no-sandbox', '--enable-unsafe-swiftshader']});
    try {
        const errors = [];
        const page = async () => {
            const p = await browser.newPage();
            p.on('pageerror', error => errors.push(error.message));
            await p.setViewport({width: 1100, height: 700});
            await p.setRequestInterception(true);
            p.on('request', request => {
                if (request.url().endsWith('/embedded-test.html')) return request.respond({status: 200, contentType: 'text/html', body: `
<!doctype html><html><head><script defer src="/static/build/3Dmol-min.js"></script><script defer src="/static/build/3Dmol.ui-min.js"></script></head><body>
<script id="atoms" type="chemical/x-xyz">2\nAtoms\nC -3 0 0\nO 3 0 0\n</script>
<div id="primary" style="position:relative;width:500px;height:600px;display:inline-block" class="viewer_3Dmoljs" data-element="atoms" data-type="xyz" data-style="sphere" data-learning="hidden" data-learning-session="${room}" data-learning-server="${origin}"></div>
<div id="other" style="position:relative;width:500px;height:600px;display:inline-block" class="viewer_3Dmoljs" data-element="atoms" data-type="xyz" data-style="sphere" data-learning="hidden" data-learning-session="${room}-other" data-learning-server="${origin}"></div>
</body></html>`});
                // Optional local asset makes this check independent of CDN access.
                if (/socket\.io(?:\.min)?\.js(?:\?|$)/.test(request.url()) && process.env.SOCKET_IO_CLIENT) {
                    return request.respond({status:200,contentType:'text/javascript',body:fs.readFileSync(process.env.SOCKET_IO_CLIENT)});
                }
                request.continue();
            });
            await p.goto(origin + '/embedded-test.html');
            await p.waitForFunction(() => !!$3Dmol.viewers.primary?.ui?.learning);
            await p.evaluate(() => $3Dmol.viewers.primary.ui.learning.ready);
            return p;
        };
        const teacher = await page(), student = await page();
        assert.equal(await teacher.$eval('#primary .molui-learning-launch', el => el.hidden), true);
        await teacher.bringToFront();
        await teacher.focus('#primary');
        await teacher.keyboard.down('Control'); await teacher.keyboard.down('Alt'); await teacher.keyboard.press('KeyI');
        await teacher.keyboard.up('Alt'); await teacher.keyboard.up('Control');
        await teacher.waitForFunction(() => $3Dmol.viewers.primary.ui.learning.getRole() === 'instructor', {timeout:10000}).catch(async error => { console.error(await teacher.evaluate(() => ({role:$3Dmol.viewers.primary.ui.learning.getRole(),status:document.querySelector('#primary .molui-learning').textContent,active:document.activeElement?.id}))); throw error; });
        await student.waitForFunction(() => $3Dmol.viewers.primary.ui.learning.getRole() === 'student', {polling:100});
        assert.deepEqual(await teacher.evaluate(() => $3Dmol.viewers.primary.ui.getControls()), {edit:true, trajectory:true});
        assert.deepEqual(await student.evaluate(() => $3Dmol.viewers.primary.ui.getControls()), {edit:false, trajectory:false});
        assert.equal(await student.evaluate(() => $3Dmol.viewers.other.ui.learning.getRole()), 'idle');
        await teacher.evaluate(() => {
            const viewer = $3Dmol.viewers.primary, state = viewer.ui;
            viewer.addLabel('First residue', {position:{x:-3,y:0,z:0}});
            viewer.addLabel('Second residue', {position:{x:3,y:0,z:0}});
            state.addDistance({x:-3,y:0,z:0},{x:3,y:0,z:0},'Two atoms',{labelText:'Contact'});
            state.setViewSettings({background:'#f0e0d0'});
        });
        await student.waitForFunction(() => $3Dmol.viewers.primary.getInternalState().labels.filter(l => !l.style?._sessionTransient).length === 3, {polling:100,timeout:10000}).catch(async error=>{ for (const [name,p] of [['teacher',teacher],['student',student]]) console.error(name, await p.evaluate(()=>{const v=$3Dmol.viewers.primary,c=v.ui.learning;return{labels:v.getInternalState().labels,role:c.getRole(),dirty:c.dirty,revision:c.revision,lastSent:c.lastSent?.labels,received:c.received?.labels,error:c.element.querySelector('[role=alert]').textContent,connected:c.socket?.connected};}));throw error; });
        const late = await page();
        await late.waitForFunction(() => $3Dmol.viewers.primary.ui.learning.getRole() === 'student' && $3Dmol.viewers.primary.getInternalState().labels.length === 3, {polling:100});
        const action = async (p, text) => p.evaluate(text => [...document.querySelectorAll('#primary .molui-learning button')].find(b => b.textContent === text).click(), text);
        await action(teacher, 'Query Atoms');
        await student.waitForFunction(() => $3Dmol.viewers.primary.selectedAtoms({}).every(a => a.clickable), {polling:100});
        await student.evaluate(() => { const atom=$3Dmol.viewers.primary.selectedAtoms({})[0]; atom.callback(atom); });
        await teacher.waitForFunction(() => document.querySelector('#primary .molui-learning-stats').textContent.includes('Responses: 1'), {polling:100});
        await action(teacher, 'Show Results');
        await teacher.waitForFunction(() => $3Dmol.viewers.primary.getInternalState().labels.some(l => l.text === '1' && l.style?._sessionTransient), {polling:100});
        await action(teacher, 'Clear Labels');
        await student.waitForFunction(() => $3Dmol.viewers.primary.selectedAtoms({}).every(a => !a.clickable), {polling:100});
        await teacher.evaluate(() => $3Dmol.viewers.primary.ui.learning.leave());
        await student.waitForFunction(() => $3Dmol.viewers.primary.ui.learning.getRole() === 'idle', {polling:100});
        assert.deepEqual(await student.evaluate(() => $3Dmol.viewers.primary.ui.getControls()), {edit:false,trajectory:false});
        await teacher.evaluate(() => $3Dmol.viewers.primary.destroy());
        assert.equal(await teacher.$('#primary .molui-learning'), null);
        assert.deepEqual(errors, []);
        console.log('Embedded hidden activation, room isolation, late joins, scene labels, queries, readonly controls, and cleanup passed.');
        // The standalone viewer uses the same controller with its URL adapter.
        await late.goto(origin + '/static/viewer.html');
        await late.waitForFunction(() => !!window.viewerApp?.learning);
        await late.evaluate(() => viewerApp.ready);
        await late.click('.molui-learning-launch');
        assert.equal(await late.$eval('.molui-learning-dialog', el => el.hidden), false);
        await late.evaluate(() => viewerApp.destroy());
        assert.equal(await late.$('.molui-learning'), null);
        console.log('Standalone classroom controls and disposal passed.');
    } finally { await browser.close(); }
})().catch(error => { console.error(error); process.exitCode=1; });
