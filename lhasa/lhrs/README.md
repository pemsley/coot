# Seed Quickstart

> Basic Rust-only template for your new Seed app.

## 1. Create a new project

1. You can use [cargo generate](https://github.com/ashleygwilliams/cargo-generate) to use this template.

    ```bash
    $ cargo generate --git https://github.com/seed-rs/seed-quickstart.git --name my-project
    $ cd my-project
    ```

1. Alternatively, simply click on the green button **Use this template** on the GitHub [profile](https://github.com/seed-rs/seed-quickstart) of this quickstart.

1. Make sure Git doesn't automatically convert your newlines to CRLF because linters don't like it.
    - Run `$ git config --global core.autocrlf` in your terminal and it should return `input` or `false`. See [Git docs](https://git-scm.com/book/en/v2/Customizing-Git-Git-Configuration) for more info.

1. Clone your new repository to your local machine. I use [GitKraken](https://www.gitkraken.com/), but you are probably a better developer than me - use your favorite terminal.

## 2. Install / check required tools

1. Make sure you have basic tools installed:

   - [Rust](https://www.rust-lang.org)
     - Check: `$ rustc -V` => `rustc 1.43.1 (8d69840ab 2020-05-04)`
     - Install: https://www.rust-lang.org/tools/install
   - [cargo-make](https://sagiegurari.github.io/cargo-make/)
     - Check: `$ cargo make -V` => `cargo-make 0.30.7`
     - Install: `$ cargo install cargo-make`

1. Platform-specific tools like `ssl` and `pkg-config`:
    - Follow recommendations in build errors (during the next chapter).
    - _Note_: Don't hesitate to write notes or a tutorial for your platform and create a PR .

## 3. Prepare your project for work

1. Open the project in your favorite IDE (I recommend [VS Code](https://code.visualstudio.com/) + [Rust Analyzer](https://rust-analyzer.github.io/)).
1. Open a new terminal tab / window and run: `cargo make serve`
1. Open a second terminal tab and run: `cargo make watch`
1. If you see errors, try to fix them or write on our [chat](https://discord.gg/JHHcHp5) or [forum](https://seed.discourse.group/).
1. Modify files like `README.md` and `Cargo.toml` as you wish.

## 4. Write your website

1. Open [localhost:8000](http://localhost:8000) in a browser (I recommend Firefox and Chrome).
1. Modify source files (e.g. `/src/lib.rs` or `/index.html`).
1. Watch compilation in the terminal tab where you run `cargo make watch`.
1. You can watch dev-server responses in the tab where you run `cargo make serve`.
1. Refresh your browser and see changes.
1. Go to step 2.

## 5. Prepare your project for deploy

1. Run `cargo make verify` in your terminal to format and lint the code.
1. Run `cargo make build_release`.
1. Upload `index.html` and `pkg` into your server's public folder.
   - Don't forget to upload also configuration files for your hosting, see the [Netlify](https://www.netlify.com/) one below.

```toml
# netlify.toml
[[redirects]]
  from = "/*"
  to = "/index.html"
  status = 200
```

## Other Seed quickstarts and projects

- [seed-rs/awesome-seed-rs](https://github.com/seed-rs/awesome-seed-rs)

---

**!!! New Rust-only quickstart in development! => [Seeder](https://github.com/MartinKavik/seeder) !!!**

---
