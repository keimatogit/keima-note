var fade = 0;

$(function() {

    $('.header-hamburger i').on('click', function() { 
        $(this).toggleClass('fa-bars');
        $(this).toggleClass('fa-times');
        $('.header-nav').fadeToggle(fade);
        $('.header').toggleClass('open');
        $('body').toggleClass('fixed'); //メニューオープン時に全体を固定
    });


    // toc(mobile)
    $('.toc-hamburger-mobile i').on('click', function() {
        $(this).toggleClass('fa-bars');
        $(this).toggleClass('fa-times');
        $('.toc-mobile').fadeToggle(fade);
        $('.toc-wrapper').toggleClass('header-open');
    });

    // toc(mobile) - メニュークリックでメニューを消す
    // ここで.header_open a など後でつけたクラスを使っても効かないのでtoc-mobileとtoc-pcを個別につくっている
    $('.toc-mobile a').click(function() { // toc内メニューをクリックした時にメニューを消す
        $('.toc-hamburger-mobile i').toggleClass('fa-bars');
        $('.toc-hamburger-mobile i').toggleClass('fa-times');
        $('.toc-mobile').fadeToggle(fade);
        $('.toc-wrapper').toggleClass('header-open');
    });

    // toc(pc)
    $('.toc-hamburger-pc i').on('click', function() {
        $(this).toggleClass('fa-bars');
        $(this).toggleClass('fa-times');
        $('.toc-pc').fadeToggle(fade);
        $('.toc-wrapper').toggleClass('sidebar-open');
        $('.toc-wrapper').toggleClass('sidebar-close');
        $('.main-wrapper').toggleClass('sidebar-open');
        $('.main-wrapper').toggleClass('sidebar-close');
    });

});


