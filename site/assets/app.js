/* Same local preferences and functional nav controls as PMx Explain / MultiOmics.
 * No analytics, third-party requests, patient data or dynamic model execution.
 */
(() => {
  const root = document.documentElement;
  const nav = document.querySelector('#site-nav');
  const burger = document.querySelector('.burger');
  const themeButton = document.querySelector('#theme-switch');
  const languages = ['en', 'fr'];
  const modes = ['system', 'light', 'dark'];
  const menus = [...document.querySelectorAll('details.goal')];
  const saved = key => { try { return localStorage.getItem(key); } catch { return null; } };
  const store = (key, value) => { try { localStorage.setItem(key, value); } catch {} };
  const changeQueryLang = lang => {
    const url = new URL(location.href);
    if (lang === 'fr') url.searchParams.delete('lang');
    else url.searchParams.set('lang', lang);
    if (url.href !== location.href) history.replaceState(history.state, '', url);
  };
  function setLanguage(value) {
    const lang = languages.includes(value) ? value : 'fr';
    root.lang = lang;
    document.querySelectorAll('[data-fr][data-en]').forEach(el => {
      el.textContent = el.getAttribute('data-' + lang);
    });
    document.querySelectorAll('[data-set-lang]').forEach(button => {
      const active = button.getAttribute('data-set-lang') === lang;
      button.classList.toggle('active', active);
      button.setAttribute('aria-pressed', String(active));
    });
    store('pharmacometrie-language', lang);
    changeQueryLang(lang);
    updateThemeLabel();
  }
  const icon = {
    system: '<rect x="3" y="4" width="18" height="13" rx="2"/><path d="M8 20h8"/>',
    light: '<circle cx="12" cy="12" r="4.2"/><path d="M12 2.5v2M12 19.5v2M2.5 12h2M19.5 12h2M5.2 5.2l1.4 1.4M17.4 17.4l1.4 1.4M18.8 5.2l-1.4 1.4M6.6 17.4l-1.4 1.4"/>',
    dark: '<path d="M20 14.5A8.5 8.5 0 1 1 9.5 4a6.8 6.8 0 0 0 10.5 10.5Z"/>'
  };
  let theme = 'system';
  function updateThemeLabel() {
    const value = root.lang === 'en'
      ? 'Theme: ' + theme : 'Thème : ' + (theme === 'system' ? 'système' : theme === 'light' ? 'clair' : 'sombre');
    themeButton?.setAttribute('aria-label', value);
    themeButton?.setAttribute('title', value);
  }
  function applyTheme(value) {
    theme = modes.includes(value) ? value : 'system';
    if (theme === 'system') root.removeAttribute('data-theme');
    else root.setAttribute('data-theme', theme);
    if (themeButton) {
      const svg = themeButton.querySelector('svg');
      if (svg) svg.innerHTML = icon[theme];
    }
    store('pk-theme', theme);
    updateThemeLabel();
  }
  function closeNavigation() {
    menus.forEach(d => { d.open = false; });
    nav?.classList.remove('open');
    burger?.setAttribute('aria-expanded','false');
  }
  const queryLang = new URLSearchParams(location.search).get('lang');
  setLanguage(languages.includes(queryLang) ? queryLang : saved('pharmacometrie-language') || 'fr');
  applyTheme(saved('pk-theme') || 'system');
  document.querySelectorAll('[data-set-lang]').forEach(button => button.addEventListener(
    'click', () => setLanguage(button.getAttribute('data-set-lang'))
  ));
  themeButton?.addEventListener('click', () => applyTheme(
    modes[(modes.indexOf(theme) + 1) % modes.length]
  ));
  burger?.addEventListener('click', () => {
    const open = nav?.classList.toggle('open');
    burger.setAttribute('aria-expanded', String(Boolean(open)));
  });
  document.addEventListener('click', event => {
    if (!event.target.closest('.site-header')) closeNavigation();
  });
  document.addEventListener('keydown', event => {
    if (event.key === 'Escape') {
      const focused = menus.find(d => d.open);
      const returnFocus = focused?.querySelector('summary') || burger;
      closeNavigation();
      returnFocus?.focus();
    }
  });
  nav?.querySelectorAll('a').forEach(link => link.addEventListener('click', closeNavigation));
})();
