/* Lightweight shared navigation, no external analytics and no patient data. */
(() => {
  const root = document.documentElement;
  const langButton = document.querySelector('#lang-switch');
  const themeButton = document.querySelector('#theme-switch');
  const burger = document.querySelector('.burger');
  const nav = document.querySelector('#site-nav');
  const menus = [...document.querySelectorAll('details.goal')];
  const languages = new Set(['fr','en']);
  const themes = new Set(['dark','light']);
  const safeGet = key => { try { return localStorage.getItem(key); } catch { return null; } };
  const safeSet = (key,value) => { try { localStorage.setItem(key,value); } catch {} };
  function setLanguage(value) {
    const lang = languages.has(value) ? value : 'fr';
    root.lang = lang;
    document.querySelectorAll('[data-fr][data-en]').forEach(el => {
      el.textContent = el.getAttribute('data-' + lang);
    });
    if(langButton) {
      langButton.setAttribute('aria-pressed',String(lang === 'en'));
      langButton.setAttribute('aria-label',lang === 'fr' ? 'Switch to English' : 'Passer en français');
    }
    safeSet('pmx-ecosystem-language',lang);
  }
  function setTheme(value) {
    const theme = themes.has(value) ? value : 'light';
    root.setAttribute('data-theme',theme);
    themeButton?.setAttribute('aria-pressed',String(theme === 'dark'));
    themeButton?.setAttribute('aria-label',theme === 'dark'?'Passer au thème clair':'Passer au thème sombre');
    safeSet('pmx-ecosystem-theme',theme);
  }
  setLanguage(safeGet('pmx-ecosystem-language') || (/^en\b/i.test(navigator.language) ? 'en':'fr'));
  setTheme(safeGet('pmx-ecosystem-theme') || (window.matchMedia?.('(prefers-color-scheme: dark)').matches ? 'dark':'light'));
  langButton?.addEventListener('click',() => setLanguage(root.lang === 'fr' ? 'en':'fr'));
  themeButton?.addEventListener('click',() => setTheme(root.dataset.theme === 'dark' ? 'light':'dark'));
  const close = () => {
    menus.forEach(d => d.open = false);
    nav?.classList.remove('open');
    burger?.setAttribute('aria-expanded','false');
  };
  burger?.addEventListener('click',() => {
    const opened = nav?.classList.toggle('open');
    burger.setAttribute('aria-expanded',String(Boolean(opened)));
  });
  document.addEventListener('click',ev => {
    if(!ev.target.closest('.site-header')) close();
  });
  document.addEventListener('keydown',ev => {
    if(ev.key === 'Escape'){
      const focused = menus.find(d => d.open);
      const backTo = focused?.querySelector('summary') || burger;
      close();
      backTo?.focus();
    }
  });
  nav?.querySelectorAll('a').forEach(link=>link.addEventListener('click',close));
})();
